/** @file pyobject_ptr.h
 *
 *  Small RAII helper for owned Python references used by the Pynac bridge.
 */

#ifndef __GINAC_PYOBJECT_PTR_H__
#define __GINAC_PYOBJECT_PTR_H__

#include <Python.h>

#include <memory>

namespace GiNaC {
namespace internal {

struct pyobject_deleter
{
	void operator()(PyObject* object) const noexcept
	{
		Py_XDECREF(object);
	}
};

using pyobject_ptr = std::unique_ptr<PyObject, pyobject_deleter>;

} // namespace internal
} // namespace GiNaC

#endif // __GINAC_PYOBJECT_PTR_H__
