get_filename_component(SUB_DIRECTORY "${CMAKE_CURRENT_LIST_DIR}" NAME)

set(TEST_FILES ElwinModelTest.h ElwinPresenterTest.h IqtModelTest.h MomentsModelTest.h SqwModelTest.h
               SymmetriseModelTest.h
)

list(TRANSFORM TEST_FILES PREPEND ${SUB_DIRECTORY}/)

set(ALL_TEST_FILES
    ${ALL_TEST_FILES} ${TEST_FILES}
    PARENT_SCOPE
)
