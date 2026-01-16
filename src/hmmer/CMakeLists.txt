# --- Handle config.h generation *a la* ./configure ----------------------------

include(CheckIncludeFile)
include(CheckFunctionExists)

CHECK_INCLUDE_FILE("netinet/in.h" HAVE_NETINET_IN_H)
CHECK_INCLUDE_FILE("sys/param.h" HAVE_SYS_PARAM_H)
CHECK_INCLUDE_FILE("sys/sysctl.h" HAVE_SYS_SYSCTL_H)

set(HMMER_DATE "Aug 2023")
set(HMMER_COPYRIGHT "Copyright (C) 2023 Howard Hughes Medical Institute.")
set(HMMER_LICENSE  "Freely distributed under the BSD open source license.")
set(HMMER_VERSION "3.4")
set(HMMER_URL "http://hmmer.org/")

configure_file(p7_config.h.cmake p7_config.h)

# --- Manually declare required sources only -----------------------------------

set(HMMER_SOURCE_DIR
    "${PROJECT_SOURCE_DIR}/vendor/hmmer/src"
)

set(HMMER_SOURCES
    build.c 
    cachedb.c 
    cachedb_shard.c 
    emit.c 
    errors.c 
    evalues.c 
    eweight.c 
    generic_decoding.c 
    generic_fwdback.c 
    generic_fwdback_chk.c 
    generic_fwdback_banded.c 
    generic_null2.c 
    generic_msv.c 
    generic_optacc.c 
    generic_stotrace.c 
    generic_viterbi.c 
    generic_vtrace.c 
    h2_io.c 
    heatmap.c 
    hmmdmstr.c 
    hmmdmstr_shard.c 
    hmmd_search_status.c 
    hmmdwrkr.c 
    hmmdwrkr_shard.c 
    hmmdutils.c 
    hmmer.c 
    logsum.c 
    modelconfig.c 
    modelstats.c 
    mpisupport.c 
    seqmodel.c 
    tracealign.c 
    p7_alidisplay.c 
    p7_bg.c 
    p7_builder.c 
    p7_domain.c 
    p7_domaindef.c 
    p7_gbands.c 
    p7_gmx.c 
    p7_gmxb.c 
    p7_gmxchk.c 
    p7_hit.c 
    p7_hmm.c 
    p7_hmmcache.c 
    p7_hmmd_search_stats.c 
    p7_hmmfile.c 
    p7_hmmwindow.c 
    p7_pipeline.c 
    p7_prior.c 
    p7_profile.c 
    p7_spensemble.c 
    p7_tophits.c 
    p7_trace.c 
    p7_scoredata.c 
    hmmpgmd2msa.c 
    fm_alphabet.c 
    fm_general.c 
    fm_sse.c 
    fm_ssv.c
)

set(HMMER_HEADERS
    cachedb.h
    cachedb_shard.h
    hmmer.h
    hmmpgmd.h
    hmmpgmd_shard.h
    p7_gbands.h
    p7_gmxb.h
    p7_gmxchk.h
    p7_hmmcache.h
)

set(HMMER_PATCHED_SOURCES
    ${CMAKE_CURRENT_BINARY_DIR}/p7_config.h
)

if(HMMER_IMPL STREQUAL "SSE")
    set(HMMER_SOURCES
        ${HMMER_SOURCES}
        impl_sse/decoding.c
        impl_sse/fwdback.c
        impl_sse/io.c
        impl_sse/mpi.c
        impl_sse/msvfilter.c
        impl_sse/null2.c
        impl_sse/optacc.c
        impl_sse/p7_omx.c
        impl_sse/p7_oprofile.c
        impl_sse/ssvfilter.c
        impl_sse/stotrace.c
        impl_sse/vitfilter.c
    )
    set(HMMER_HEADERS
        ${HMMER_HEADERS}
        impl_sse/impl_sse.h
    )
elseif(HMMER_IMPL STREQUAL "NEON")
    set(HMMER_SOURCES
        ${HMMER_SOURCES}
        impl_neon/decoding.c
        impl_neon/fwdback.c
        impl_neon/io.c
        impl_neon/mpi.c
        impl_neon/msvfilter.c
        impl_neon/null2.c
        impl_neon/optacc.c
        impl_neon/p7_omx.c
        impl_neon/p7_oprofile.c
        impl_neon/ssvfilter.c
        impl_neon/stotrace.c
        impl_neon/vitfilter.c
    )
    set(HMMER_HEADERS
        ${HMMER_HEADERS}
        impl_neon/impl_neon.h
    )
elseif(HMMER_IMPL STREQUAL "VMX")
    set(HMMER_SOURCES
        ${HMMER_SOURCES}
        impl_vmx/decoding.c
        impl_vmx/fwdback.c
        impl_vmx/io.c
        impl_vmx/mpi.c
        impl_vmx/msvfilter.c
        impl_vmx/null2.c
        impl_vmx/optacc.c
        impl_vmx/p7_omx.c
        impl_vmx/p7_oprofile.c
        impl_vmx/stotrace.c
        impl_vmx/vitfilter.c
    )
    set(HMMER_HEADERS
        ${HMMER_HEADERS}
        impl_vmx/impl_vmx.h
    )
else()
    message(FATAL_ERROR "No HMMER implementation defined")
endif()

# --- Copy/patch source files if needed ----------------------------------------

foreach(_file IN ITEMS ${HMMER_SOURCES} ${HMMER_HEADERS})
    if(EXISTS ${PROJECT_SOURCE_DIR}/patches/${_file}.patch)
        if(NOT EXISTS ${HMMER_SOURCE_DIR}/${_file})
            add_custom_command(
                OUTPUT
                    ${_file}
                COMMENT
                    "Creating ${_file}"
                COMMAND
                    cmake -E touch ${CMAKE_CURRENT_BINARY_DIR}/empty.txt
                COMMAND
                    ${Python_EXECUTABLE} ${PROJECT_SOURCE_DIR}/src/scripts/apply_patch.py
                        --input ${CMAKE_CURRENT_BINARY_DIR}/empty.txt
                        --patch ${PROJECT_SOURCE_DIR}/patches/${_file}.patch
                        --output ${CMAKE_CURRENT_BINARY_DIR}/${_file}
                DEPENDS
                    ${PROJECT_SOURCE_DIR}/patches/${_file}.patch
            )
        else()
            add_custom_command(
                OUTPUT
                    ${_file}
                COMMENT
                    "Patching ${_file}"
                COMMAND
                    ${Python_EXECUTABLE} ${PROJECT_SOURCE_DIR}/src/scripts/apply_patch.py
                        --input ${HMMER_SOURCE_DIR}/${_file}
                        --patch ${PROJECT_SOURCE_DIR}/patches/${_file}.patch
                        --output ${CMAKE_CURRENT_BINARY_DIR}/${_file}
                DEPENDS
                    ${HMMER_SOURCE_DIR}/${_file}
                    ${PROJECT_SOURCE_DIR}/patches/${_file}.patch
            )
        endif()
    else()
        add_custom_command(
            OUTPUT
                ${_file}
            COMMENT
                "Copying ${_file}"
            COMMAND
                cmake -E copy ${HMMER_SOURCE_DIR}/${_file} ${CMAKE_CURRENT_BINARY_DIR}/${_file}
            DEPENDS
                ${HMMER_SOURCE_DIR}/${_file}
        )
    endif()
    set(HMMER_PATCHED_SOURCES ${HMMER_PATCHED_SOURCES} ${CMAKE_CURRENT_BINARY_DIR}/${_file})
endforeach()

# --- Build libhmmer with copied/patched sources -------------------------------

add_library(libhmmer ${HMMER_SOURCES} ${HMMER_HEADERS})
target_link_libraries(libhmmer PUBLIC libeasel)
target_include_directories(libhmmer PUBLIC ${CMAKE_CURRENT_BINARY_DIR})
target_include_directories(libhmmer PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})

# --- Install shared libraries & headers if required ---------------------------

if(DEFINED PYHMMER_INSTALL_LIBS_DIR)
    install(TARGETS libhmmer DESTINATION "${SKBUILD_PLATLIB_DIR}/${PYHMMER_INSTALL_LIBS_DIR}")
    if(CMAKE_SYSTEM_NAME STREQUAL "Darwin")
        set_target_properties(libhmmer PROPERTIES INSTALL_RPATH "@loader_path/")
    else()
        set_target_properties(libhmmer PROPERTIES INSTALL_RPATH "\$ORIGIN/")
    endif()    
    cmake_path(APPEND SKBUILD_PLATLIB_DIR "${PYHMMER_INSTALL_LIBS_DIR}" "include" "libhmmer" OUTPUT_VARIABLE INCLUDE_DEST_FOLDER)
    foreach(_file IN ITEMS ${HMMER_PATCHED_SOURCES})
        cmake_path(GET _file EXTENSION LAST_ONLY _ext)
        if(_ext STREQUAL ".h")
            cmake_path(RELATIVE_PATH _file BASE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}" OUTPUT_VARIABLE _local)
            cmake_path(GET _local PARENT_PATH _dest)
            install(FILES "${_file}" DESTINATION "${INCLUDE_DEST_FOLDER}/${_dest}")
        endif()
    endforeach()
endif()
