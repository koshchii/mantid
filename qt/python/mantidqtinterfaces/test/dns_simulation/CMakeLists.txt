# Tests for DNS GUI

set(TEST_PY_FILES
    simulation_sub_common_presenter_test.py
    simulation_sub_powder_widget_test.py
    simulation_sub_powder_presenter_test.py
    simulation_sub_sc_model_test.py
    simulation_sub_table_model_test.py
    simulation_widget_test.py
    simulation_helpers_test.py
    simulation_model_test.py
    simulation_presenter_test.py
    simulation_sub_powder_model_test.py
    simulation_sub_table_widget_test.py
    simulation_sub_sc_presenter_test.py
    simulation_sub_sc_widget_test.py
    simulation_sub_table_presenter_test.py
)

check_tests_valid(${CMAKE_CURRENT_SOURCE_DIR} ${TEST_PY_FILES})

# Prefix for test name=PythonAlgorithms
set(PYUNITTEST_QT_API pyqt5) # force to use qt5
pyunittest_add_test(${CMAKE_CURRENT_SOURCE_DIR} python.DNSReductionQt5 ${TEST_PY_FILES})
unset(PYUNITTEST_QT_API)
