if(NOT DEFINED PROGRAM)
  message(FATAL_ERROR "PROGRAM not set")
endif()

if(NOT DEFINED WORKING_DIR)
  message(FATAL_ERROR "WORKING_DIR not set")
endif()

if(NOT DEFINED SCRIPT)
  message(FATAL_ERROR "SCRIPT not set")
endif()

if(NOT DEFINED OUTPUT_FILE)
  message(FATAL_ERROR "OUTPUT_FILE not set")
endif()

if(NOT DEFINED SUMMARY_COMPONENT)
  message(FATAL_ERROR "SUMMARY_COMPONENT not set")
endif()

if(NOT DEFINED EVOLUTION_COMPONENT)
  message(FATAL_ERROR "EVOLUTION_COMPONENT not set")
endif()

set(output_path "${WORKING_DIR}/${OUTPUT_FILE}")
file(REMOVE "${output_path}")

execute_process(
  COMMAND ${CMAKE_COMMAND} -E env OMP_NUM_THREADS=1 ${PROGRAM} ${SCRIPT}
  WORKING_DIRECTORY "${WORKING_DIR}"
  RESULT_VARIABLE result
)

if(NOT result EQUAL 0)
  message(FATAL_ERROR "Command failed with exit code ${result}")
endif()

if(NOT EXISTS "${output_path}")
  message(FATAL_ERROR "Expected output file not found: ${output_path}")
endif()

file(STRINGS "${output_path}" content_lines)

set(found_summary FALSE)
set(found_evolution FALSE)

foreach(line IN LISTS content_lines)
  if(line MATCHES "^summary,")
    string(FIND "${line}" "${SUMMARY_COMPONENT}" summary_pos)
    if(NOT summary_pos EQUAL -1)
      set(found_summary TRUE)
    endif()
  endif()

  if(line MATCHES "^evolution,")
    string(FIND "${line}" "${EVOLUTION_COMPONENT}" evolution_pos)
    if(NOT evolution_pos EQUAL -1)
      set(found_evolution TRUE)
    endif()
  endif()
endforeach()

if(NOT found_summary)
  message(FATAL_ERROR "Missing summary DIB rows in ${output_path}")
endif()

if(NOT found_evolution)
  message(FATAL_ERROR "Missing evolution DIB rows in ${output_path}")
endif()
