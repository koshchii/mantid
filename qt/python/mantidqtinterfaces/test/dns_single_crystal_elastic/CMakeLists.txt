# Tests for DNS GUI

set(TEST_PY_FILES
    data_structures/dns_dialog_test.py
    data_structures/dns_single_crystal_map_test.py
    options/elastic_single_crystal_options_model_test.py
    options/elastic_single_crystal_options_presenter_test.py
    options/elastic_single_crystal_options_widget_test.py
    plot/elastic_single_crystal_plot_model_test.py
    plot/elastic_single_crystal_plot_widget_test.py
    plot/elastic_single_crystal_plot_presenter_test.py
    script_generator/elastic_single_crystal_script_generator_presenter_test.py
    script_generator/elastic_single_crystal_script_generator_widget_test.py
    script_generator/elastic_single_crystal_script_generator_model_test.py
    scripts/md_single_crystal_elastic_test.py
)

check_tests_valid(${CMAKE_CURRENT_SOURCE_DIR} ${TEST_PY_FILES})

# Prefix for test name=PythonAlgorithms
set(PYUNITTEST_QT_API pyqt5) # force to use qt5
pyunittest_add_test(${CMAKE_CURRENT_SOURCE_DIR} python.DNSReductionQt5 ${TEST_PY_FILES})
unset(PYUNITTEST_QT_API)
