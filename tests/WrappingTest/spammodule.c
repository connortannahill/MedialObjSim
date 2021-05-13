#define PY_SSIZE_T_CLEAN
#define NULL 0
#include <Python/Python.h>

// Defining the new exception that is unique to this module
static PyObject *SpamError;

static PyObject * spam_system(PyObject *self, PyObject *args) {
    const char *command;
    int sts;

    if (!PyArg_ParseTuple(args, "s", &command))
        return NULL;

    sts = system(command);
    return PyLong_FromLong(sts);
}

PyMODINIT_FUNC PyInit_spam(void) {
    PyObject *m;

    m = PyModule_Create(&spammodule);

    if (m == NULL) 
        return NULL;

    SpamError = PyErr_NewException("spam.error", NULL, NULL);
    Py_INCREF(SpamError);
    PyModule_AddObject(m, "error", SpamError);
    return m;
}
