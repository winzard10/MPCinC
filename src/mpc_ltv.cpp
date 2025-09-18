#include "mpc_ltv.hpp"
#include <OsqpEigen/OsqpEigen.h>
#include <cmath>
#include <cassert>

using Eigen::MatrixXd;
using Eigen::VectorXd;
using Eigen::SparseMatrix;
using Eigen::Triplet;

#ifndef OSQP_INFTY
#define OSQP_INFTY 1e20
#endif

static inline double clamp(double x, double a, double b){ return std::min(std::max(x,a),b); }
void LTV_MPC::angleWrap(double& a){ while(a>M_PI)a-=2*M_PI; while(a<=-M_PI)a+=2*M_PI; }

LTV_MPC::LTV_MPC(const MPCParams& p) : P(p) {
    nom_.x.resize(P.N+1);
    nom_.u.resize(P.N, Eigen::Vector2d::Zero());
}

void LTV_MPC::resetWarmStart(){
    for (auto& u : nom_.u) u.setZero();
}

void LTV_MPC::buildLinearization(const MPCRef& ref) {
    for (int k=0; k<=P.N; ++k) {
        nom_.x[k].ey    = 0.0;
        nom_.x[k].epsi  = 0.0;
        nom_.x[k].v     = (k < P.N) ? ref.hp[k].v_ref : ref.hp[P.N-1].v_ref;
        nom_.x[k].delta = 0.0;
    }
    for (int k=0; k<P.N; ++k) nom_.u[k].setZero();
}

MPCControl LTV_MPC::solve(const MPCState& x0, const MPCRef& ref){
    buildLinearization(ref);
    return solveQP(x0, ref);
}

MPCControl LTV_MPC::solveQP(const MPCState& x0, const MPCRef& ref){
    const int nx = 4;             // [ey, epsi, v, delta]
    const int nu = 2;             // [a, ddelta]
    const int N  = P.N;

    // Decision vector z = [x0..xN, u0..uN-1]
    const int NX = (N+1)*nx;
    const int NU = N*nu;
    const int NS = 2*(N+1);
    const int NZ = NX + NU + NS;

    // ---------- Cost ----------
    std::vector<Triplet<double>> Ht;
    VectorXd g = VectorXd::Zero(NZ);

    auto idx_x = [&](int k,int i){ return k*nx + i; };       // state
    auto idx_u = [&](int k,int j){ return NX + k*nu + j; };  // input

    auto idx_sig_ey    = [&](int k){ return NX + NU + (2*k + 0); };
    auto idx_sig_delta = [&](int k){ return NX + NU + (2*k + 1); };

    for (int k=0;k<N;++k){
        // state tracking
        Ht.emplace_back(idx_x(k,0), idx_x(k,0), P.wy);
        Ht.emplace_back(idx_x(k,1), idx_x(k,1), P.wpsi);
        Ht.emplace_back(idx_x(k,2), idx_x(k,2), P.wv);
        Ht.emplace_back(idx_sig_ey(k),    idx_sig_ey(k),    P.w_sigma_ey);
        Ht.emplace_back(idx_sig_delta(k), idx_sig_delta(k), P.w_sigma_delta);
        g(idx_x(k,2)) += -2.0 * P.wv * ref.hp[k].v_ref; // linear term for v_ref

        // input effort
        Ht.emplace_back(idx_u(k,0), idx_u(k,0), P.wa);
        Ht.emplace_back(idx_u(k,1), idx_u(k,1), P.wdd);

        // input slew
        if (k>0){
            for (int j=0;j<nu;++j){
                const int uk   = idx_u(k,j);
                const int ukm1 = idx_u(k-1,j);
                const double w = (j==0? P.wda : P.wddd);
                Ht.emplace_back(uk,    uk,    w);
                Ht.emplace_back(ukm1,  ukm1,  w);
                Ht.emplace_back(uk,    ukm1, -w);
                Ht.emplace_back(ukm1,  uk,   -w);
            }
        }
    }
    // terminal weights
    Ht.emplace_back(idx_x(N,0), idx_x(N,0), P.wyf);
    Ht.emplace_back(idx_x(N,1), idx_x(N,1), P.wpsif);

    SparseMatrix<double> H(NZ, NZ);
    H.setFromTriplets(Ht.begin(), Ht.end());

    // ---------- Dynamics equalities: Aeq z = beq ----------
    const int meq = N*nx + nx; // N transitions + init
    std::vector<Triplet<double>> Aeqt;
    VectorXd beq = VectorXd::Zero(meq);

    auto row_dyn = [&](int k,int i){ return k*nx + i; };
    auto row_ic  = [&](int i){ return N*nx + i; };

    // initial condition
    beq(row_ic(0)) = x0.ey;
    beq(row_ic(1)) = x0.epsi;
    beq(row_ic(2)) = x0.v;
    beq(row_ic(3)) = x0.delta;
    for (int i=0;i<nx;++i) Aeqt.emplace_back(row_ic(i), idx_x(0,i), 1.0);

    // linearized discrete dynamics
    for (int k=0;k<N;++k){
        const double eybar   = nom_.x[k].ey;
        const double epsibar = nom_.x[k].epsi;
        const double vbar    = nom_.x[k].v;
        const double delbar  = nom_.x[k].delta;
        const double kap     = ref.hp[k].kappa;

        const double cpsi = std::cos(epsibar);
        const double spsi = std::sin(epsibar);
        const double den  = 1.0 - kap*eybar;
        const double den2 = den*den;
        const double cosd = std::cos(delbar);
        const double sec2 = 1.0 / (cosd*cosd);

        MatrixXd A = MatrixXd::Zero(nx, nx);
        MatrixXd B = MatrixXd::Zero(nx, nu);
        VectorXd d = VectorXd::Zero(nx);

        // ey_dot
        A(0,0) = 0.0;
        A(0,1) = vbar * cpsi;
        A(0,2) = spsi;
        // epsi_dot
        A(1,0) = - vbar * kap*kap / den2;
        A(1,1) = 0.0;
        A(1,2) = (1.0/P.L)*std::tan(delbar) - (kap/den);
        A(1,3) = (vbar/P.L) * sec2;
        // inputs
        B(2,0) = 1.0; // a
        B(3,1) = 1.0; // ddelta

        const Eigen::Vector4d xbar(eybar, epsibar, vbar, delbar);
        const double ey_dot    = vbar * spsi;
        const double epsi_dot  = (vbar/P.L)*std::tan(delbar) - vbar*kap/den;
        d(0) = ey_dot   - (A.row(0) * xbar)(0);
        d(1) = epsi_dot - (A.row(1) * xbar)(0);
        d(2) = 0.0;
        d(3) = 0.0;

        MatrixXd Ad = MatrixXd::Identity(nx,nx) + P.dt * A;
        MatrixXd Bd = P.dt * B;
        VectorXd cd = P.dt * d;

        for (int i=0;i<nx;++i){
            // x_{k+1}
            Aeqt.emplace_back(row_dyn(k,i), idx_x(k+1,i), 1.0);
            // -Ad x_k
            for (int j=0;j<nx;++j){
                const double val = -Ad(i,j);
                if (val!=0.0) Aeqt.emplace_back(row_dyn(k,i), idx_x(k,j), val);
            }
            // -Bd u_k
            for (int j=0;j<nu;++j){
                const double val = -Bd(i,j);
                if (val!=0.0) Aeqt.emplace_back(row_dyn(k,i), idx_u(k,j), val);
            }
            // RHS
            beq(row_dyn(k,i)) = -cd(i);
        }
    }

    SparseMatrix<double> Aeq(meq, NZ);
    Aeq.setFromTriplets(Aeqt.begin(), Aeqt.end());

    // ---------- Variable bounds as constraints ----------
    VectorXd lb = VectorXd::Constant(NZ, -OSQP_INFTY);
    VectorXd ub = VectorXd::Constant(NZ,  OSQP_INFTY);

    // inputs
    for (int k=0;k<N;++k){
        lb(idx_u(k,0)) = P.a_min;
        ub(idx_u(k,0)) = P.a_max;
        lb(idx_u(k,1)) = -P.ddelta_max;
        ub(idx_u(k,1)) =  P.ddelta_max;
    }
    // speed (hard bounds)
    for (int k=0;k<=N;++k){
        lb(idx_x(k,2)) = P.v_min;
        ub(idx_x(k,2)) = P.v_max;
    }

    // ey, delta: leave infinite (soft constraints below)

    // slacks must be nonnegative
    for (int k=0;k<=N;++k){
        lb[idx_sig_ey(k)]    = 0.0;  // 0 <= sigma
        lb[idx_sig_delta(k)] = 0.0;
        // ub already +inf

    // ---------- Build combined constraint matrix A and bounds ----------
    // Rows = equalities (meq) + soft ey (2*(N+1)) + soft delta (2*(N+1)) + identity for variable bounds (NZ)
    const int m_soft_ey    = 2*(N+1);
    const int m_soft_delta = 2*(N+1);
    const int m = meq + m_soft_ey + m_soft_delta + NZ;

    std::vector<Triplet<double>> At;
    At.reserve(Aeq.nonZeros() + m_soft_ey*3 + m_soft_delta*3 + NZ);

    // 1) copy Aeq (equalities)
    for (int col=0; col<Aeq.outerSize(); ++col){
        for (SparseMatrix<double>::InnerIterator it(Aeq,col); it; ++it){
            At.emplace_back(it.row(), it.col(), it.value());
        }
    }

    // 2) soft inequalities for ey and delta
    int row = meq;

    // ey_k - sigma_ey_k <= ey_max
    for (int k=0; k<=N; ++k, ++row){
        At.emplace_back(row, idx_x(k,0),    1.0);  // + ey_k
        At.emplace_back(row, idx_sig_ey(k), -1.0); // - sigma_ey_k
    }
    // -ey_k - sigma_ey_k <= ey_max
    for (int k=0; k<=N; ++k, ++row){
        At.emplace_back(row, idx_x(k,0),   -1.0);  // - ey_k
        At.emplace_back(row, idx_sig_ey(k), -1.0); // - sigma_ey_k
    }

    // delta_k - sigma_delta_k <= delta_max
    for (int k=0; k<=N; ++k, ++row){
        At.emplace_back(row, idx_x(k,3),        1.0);  // + delta_k
        At.emplace_back(row, idx_sig_delta(k), -1.0);  // - sigma_delta_k
    }
    // -delta_k - sigma_delta_k <= delta_max
    for (int k=0; k<=N; ++k, ++row){
        At.emplace_back(row, idx_x(k,3),       -1.0);  // - delta_k
        At.emplace_back(row, idx_sig_delta(k), -1.0);  // - sigma_delta_k
    }

    // 3) identity rows to impose variable simple bounds (lb <= z <= ub)
    for (int i=0; i<NZ; ++i, ++row){
        At.emplace_back(row, i, 1.0);
    }
    assert(row == m);

    SparseMatrix<double> A(m, NZ);
    A.setFromTriplets(At.begin(), At.end());

    // Bounds vectors
    VectorXd l(m), u(m);

    // equalities: l = u = beq
    l.head(meq) = beq;
    u.head(meq) = beq;

    // soft ey rows (upper = ey_max, lower = -inf)
    int offs = meq;
    for (int k=0; k<=N; ++k) { l(offs+k) = -OSQP_INFTY; u(offs+k) =  P.ey_max; }
    offs += (N+1);
    for (int k=0; k<=N; ++k) { l(offs+k) = -OSQP_INFTY; u(offs+k) =  P.ey_max; }
    offs += (N+1);

    // soft delta rows (upper = delta_max, lower = -inf)
    for (int k=0; k<=N; ++k) { l(offs+k) = -OSQP_INFTY; u(offs+k) =  P.delta_max; }
    offs += (N+1);
    for (int k=0; k<=N; ++k) { l(offs+k) = -OSQP_INFTY; u(offs+k) =  P.delta_max; }
    offs += (N+1);

    // identity rows for variable bounds (lb <= z <= ub)
    l.tail(NZ) = lb;
    u.tail(NZ) = ub;

    // ---------- Solve with OSQP ----------
    OsqpEigen::Solver solver;
    solver.settings()->setWarmStart(true);
    solver.settings()->setVerbosity(false);
    solver.settings()->setAbsoluteTolerance(1e-4);
    solver.settings()->setRelativeTolerance(1e-4);

    solver.data()->setNumberOfVariables(NZ);
    solver.data()->setNumberOfConstraints(m);
    solver.data()->setHessianMatrix(H);
    solver.data()->setGradient(g);
    solver.data()->setLinearConstraintsMatrix(A);
    solver.data()->setLowerBound(l);
    solver.data()->setUpperBound(u);

    if (!solver.initSolver()) return {};

    auto status = solver.solveProblem();
    if (status != OsqpEigen::ErrorExitFlag::NoError) return {};

        Eigen::VectorXd z = solver.getSolution();
        double sum_sigma_ey = 0.0, sum_sigma_delta = 0.0;
        for (int k=0; k<=P.N; ++k) {
            sum_sigma_ey    += z[idx_sig_ey(k)];
            sum_sigma_delta += z[idx_sig_delta(k)];
        }
        MPCControl out{};
        out.a      = z(idx_u(0,0));
        out.ddelta = z(idx_u(0,1));
        out.ok     = true;
        return out;
    }
}
