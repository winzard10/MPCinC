cmake_minimum_required(VERSION 3.16)
project(hwy_mpc LANGUAGES CXX)

# ===== C++ setup =====
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

if (NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
endif()

if (MSVC)
  add_compile_options(/W4 /permissive-)
else()
  add_compile_options(-Wall -Wextra -Wpedantic)
endif()

# ===== Tell CMake where you installed packages (~/.local) =====
# This is where `cmake --install` put osqp + OsqpEigen from Option A.
list(APPEND CMAKE_PREFIX_PATH "$ENV{HOME}/.local")

# ===== Dependencies (use installed packages) =====
# Eigen3
find_package(Eigen3 QUIET NO_MODULE)
if (NOT Eigen3_FOUND)
  # Fallback include path used on Ubuntu if libeigen3-dev is present
  message(STATUS "Eigen3 not found via package config; falling back to /usr/include/eigen3")
  add_library(Eigen3::Eigen INTERFACE IMPORTED)
  set_target_properties(Eigen3::Eigen PROPERTIES
    INTERFACE_INCLUDE_DIRECTORIES "/usr/include/eigen3")
endif()

# OsqpEigen (this also finds the installed osqp via its config files)
find_package(OsqpEigen REQUIRED)         # from ~/.local
# -> Provides target: OsqpEigen::OsqpEigen

# ===== Library: centerline =====
add_library(centerline STATIC
  src/centerline_map.cpp
)
target_include_directories(centerline PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include)
target_link_libraries(centerline PUBLIC Eigen3::Eigen) 

# ===== Executable: centerline_demo =====
add_executable(centerline_demo
  src/demo_centerline.cpp
)
target_link_libraries(centerline_demo PRIVATE centerline)

# ===== Executable: sim_lane_follow =====
add_executable(sim_lane_follow
  src/sim_lane_follow.cpp
)
target_link_libraries(sim_lane_follow PRIVATE centerline)

# ===== Library: mpc (LTV-MPC using OsqpEigen) =====
add_library(mpc STATIC
  src/mpc_ltv.cpp
)
target_include_directories(mpc PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include)
# NOTE: PUBLIC so Eigen’s include path propagates to sim_lane_follow (which includes mpc headers)
target_link_libraries(mpc PUBLIC Eigen3::Eigen OsqpEigen::OsqpEigen)

# Attach MPC to simulator
target_link_libraries(sim_lane_follow PRIVATE mpc)

# ===== Library: obstacles =====
add_library(obstacles STATIC
  src/obstacles.cpp
)
target_include_directories(obstacles PUBLIC
  ${CMAKE_CURRENT_SOURCE_DIR}/include
)

# Link it into the sim
target_link_libraries(sim_lane_follow PRIVATE obstacles)

# ===== Library: corridor_planner =====
add_library(corridor_planner STATIC
  src/corridor_planner.cpp
)

# Make the header(s) visible to anything that links against it
target_include_directories(corridor_planner PUBLIC
  ${CMAKE_CURRENT_SOURCE_DIR}/include
)

# Corridor planner depends on Eigen
target_link_libraries(corridor_planner PUBLIC Eigen3::Eigen)

# Link it into the sim
target_link_libraries(sim_lane_follow PRIVATE corridor_planner)

# ===== Data handling =====
set(HWY_MPC_DATA_DIR "${CMAKE_CURRENT_SOURCE_DIR}/data")

add_custom_command(TARGET centerline_demo POST_BUILD
  COMMAND ${CMAKE_COMMAND} -E make_directory "$<TARGET_FILE_DIR:centerline_demo>/data"
  COMMAND ${CMAKE_COMMAND} -E copy_directory
          "${HWY_MPC_DATA_DIR}"
          "$<TARGET_FILE_DIR:centerline_demo>/data"
  COMMENT "Copying data/ next to centerline_demo"
)

add_custom_command(TARGET sim_lane_follow POST_BUILD
  COMMAND ${CMAKE_COMMAND} -E make_directory "$<TARGET_FILE_DIR:sim_lane_follow>/data"
  COMMAND ${CMAKE_COMMAND} -E copy_directory
          "${HWY_MPC_DATA_DIR}"
          "$<TARGET_FILE_DIR:sim_lane_follow>/data"
  COMMENT "Copying data/ next to sim_lane_follow"
)

# Optional runtime macro for installed layout
target_compile_definitions(centerline_demo PRIVATE HWY_MPC_INSTALLED_DATA_SUBDIR="share/hwy_mpc/data")
target_compile_definitions(sim_lane_follow  PRIVATE HWY_MPC_INSTALLED_DATA_SUBDIR="share/hwy_mpc/data")

# ===== Install (optional) =====
install(TARGETS centerline_demo sim_lane_follow RUNTIME DESTINATION bin)
install(DIRECTORY "${HWY_MPC_DATA_DIR}/" DESTINATION share/hwy_mpc/data)

message(STATUS "Configured hwy_mpc at ${CMAKE_SOURCE_DIR}")
