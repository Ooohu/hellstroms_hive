include_directories (${CMAKE_SOURCE_DIR}/inc/)
include_directories ( ${CMAKE_SOURCE_DIR}/tinyxml/)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/xgboost/include/)

add_library (hiveLib bdt_eff.cxx bdt_datamc.cxx bdt_varplot.cxx bdt_spec.cxx bdt_sig.cxx bdt_recomc.cxx bdt_response.cxx bdt_app.cxx bdt_train.cxx bdt_file.cxx load_mva_param.cxx bdt_precalc.cxx TruncMean.cxx bdt_boxcut.cxx bdt_test.cxx bdt_vertex_eff.cxx bdt_scatter.cxx bdt_systematics.cxx bdt_var.cxx)
target_link_libraries (hiveLib tinyxmllib xgboost ${ROOT_LIBRARIES})

message(STATUS "ROOT_LIBRARIES @: ${ROOT_LIBRARIES}")


add_executable (hive hive.cxx)
target_link_libraries (hive tinyxmllib hiveLib xgboost ${ROOT_LIBRARIES})

#add_executable (validate_hive validate_hive.cxx)
#target_link_libraries (validate_hive tinyxmllib hiveLib xgboost ${ROOT_LIBRARIES})

#add_executable (validate_trk_2 validates/validate_trk_2.cxx)
#target_link_libraries (validate_trk_2 tinyxmllib hiveLib xgboost ${ROOT_LIBRARIES})


set_output_directory(hive ${CMAKE_BINARY_DIR}/bin)
#set_output_directory(validate_hive ${CMAKE_SOURCE_DIR}/build/bin)
#set_output_directory(hive ${CMAKE_CURRENT_SOURCE_DIR}/bin)
#set_output_directory(validate_hive ${CMAKE_CURRENT_SOURCE_DIR}/bin)


#Comment these in to build some lagacy code
#add_subdirectory(obsol)

#this is where quite validation scripts can live
#add_subdirectory(validates)


