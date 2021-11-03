cmake_minimum_required( VERSION 3.14 )
project(ex01)
find_package( benchmark REQUIRED )
find_package( Threads REQUIRED )

function( configure_exercise_target targetname )
	target_compile_features( ${targetname} PRIVATE cxx_std_17 )
	target_compile_options( ${targetname} PRIVATE -xSKYLAKE )
	target_link_libraries( ${targetname} PRIVATE benchmark::benchmark )
endfunction()

add_executable( matVec matVec.cpp )
configure_exercise_target( matVec )

add_executable( lambda lambda.cpp )
configure_exercise_target( lambda )

add_executable( async async.cpp )
configure_exercise_target ( async )
target_link_libraries( async PRIVATE Threads::Threads )
