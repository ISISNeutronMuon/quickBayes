set(Qldata_module_name "QLdata")
set(QLdata_main ${CMAKE_CURRENT_SOURCE_DIR}/QLdata_main.f90)
set(QLdata_only_arg "qldata")
list(APPEND QLdata_other_sources 
${CMAKE_CURRENT_SOURCE_DIR}/QLdata_subs.f90
${CMAKE_CURRENT_SOURCE_DIR}/Util.f90
${CMAKE_CURRENT_SOURCE_DIR}/Bayes.f90 
${CMAKE_CURRENT_SOURCE_DIR}/Four.f90)
f2py_fortran_library(${Qldata_module_name} ${QLdata_main} ${QLdata_only_arg} ${QLdata_other_sources})

set(QLres_module_name "QLres")
set(QLres_main ${CMAKE_CURRENT_SOURCE_DIR}/QLres_main.f90)
set(QLres_only_arg "qlres")
list(APPEND QLres_other_sources 
${CMAKE_CURRENT_SOURCE_DIR}/QLres_subs.f90
${CMAKE_CURRENT_SOURCE_DIR}/BlrRes.f90
${CMAKE_CURRENT_SOURCE_DIR}/Bayes.f90 
${CMAKE_CURRENT_SOURCE_DIR}/Four.f90 
${CMAKE_CURRENT_SOURCE_DIR}/Util.f90
)
f2py_fortran_library(${QLres_module_name} ${QLres_main} ${QLres_only_arg} ${QLres_other_sources})

set(QLse_module_name "QLse")
set(QLse_main ${CMAKE_CURRENT_SOURCE_DIR}/QLse_main.f90)
set(QLse_only_arg "qlstexp")
list(APPEND QLse_other_sources
${CMAKE_CURRENT_SOURCE_DIR}/QLse_subs.f90
${CMAKE_CURRENT_SOURCE_DIR}/BlrRes.f90
${CMAKE_CURRENT_SOURCE_DIR}/Bayes.f90
${CMAKE_CURRENT_SOURCE_DIR}/Four.f90
${CMAKE_CURRENT_SOURCE_DIR}/Util.f90
${CMAKE_CURRENT_SOURCE_DIR}/Simopt.f90
)
f2py_fortran_library(${QLse_module_name} ${QLse_main} ${QLse_only_arg} ${QLse_other_sources})

set(Quest_module_name "Quest")
set(Quest_main ${CMAKE_CURRENT_SOURCE_DIR}/Quest_main.f90)
set(Quest_only_arg "quest")
list(APPEND Quest_other_sources 
${CMAKE_CURRENT_SOURCE_DIR}/Quest_subs.f90
${CMAKE_CURRENT_SOURCE_DIR}/BlrRes.f90
${CMAKE_CURRENT_SOURCE_DIR}/Bayes.f90
${CMAKE_CURRENT_SOURCE_DIR}/Four.f90 
${CMAKE_CURRENT_SOURCE_DIR}/Util.f90
${CMAKE_CURRENT_SOURCE_DIR}/Simopt.f90
)
f2py_fortran_library(${Quest_module_name} ${Quest_main} ${Quest_only_arg} ${Quest_other_sources})

set(ResNorm_module_name "ResNorm")
set(ResNorm_main ${CMAKE_CURRENT_SOURCE_DIR}/ResNorm_main.f90)
set(ResNorm_only_arg "resnorm")
list(APPEND ResNorm_other_sources 
${CMAKE_CURRENT_SOURCE_DIR}/ResNorm_subs.f90
${CMAKE_CURRENT_SOURCE_DIR}/BlrRes.f90
${CMAKE_CURRENT_SOURCE_DIR}/Bayes.f90
${CMAKE_CURRENT_SOURCE_DIR}/Four.f90 
${CMAKE_CURRENT_SOURCE_DIR}/Util.f90
)
f2py_fortran_library(${ResNorm_module_name} ${ResNorm_main} ${ResNorm_only_arg} ${ResNorm_other_sources})