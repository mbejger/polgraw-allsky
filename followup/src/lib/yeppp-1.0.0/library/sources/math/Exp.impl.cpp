/*
 *                       Yeppp! library implementation
 *                   This file is auto-generated by Peach-Py,
 *        Portable Efficient Assembly Code-generator in Higher-level Python,
 *                  part of the Yeppp! library infrastructure
 * This file is part of Yeppp! library and licensed under the New BSD license.
 * See LICENSE.txt for the full text of the license.
 */

#include <yepBuiltin.h>
#include <yepMath.h>


extern "C" YEP_LOCAL_SYMBOL YepStatus _yepMath_Exp_V64f_V64f_Default(const Yep64f *YEP_RESTRICT xPointer, Yep64f *YEP_RESTRICT yPointer, YepSize length) {
	if YEP_UNLIKELY(xPointer == YEP_NULL_POINTER) {
		return YepStatusNullPointer;
	}
	if YEP_UNLIKELY(yepBuiltin_GetPointerMisalignment(xPointer, sizeof(Yep64f)) != 0) {
		return YepStatusMisalignedPointer;
	}
	if YEP_UNLIKELY(yPointer == YEP_NULL_POINTER) {
		return YepStatusNullPointer;
	}
	if YEP_UNLIKELY(yepBuiltin_GetPointerMisalignment(yPointer, sizeof(Yep64f)) != 0) {
		return YepStatusMisalignedPointer;
	}
	while (length-- != 0) {
		const Yep64f x = *xPointer++;
		const Yep64f y = yepBuiltin_Exp_64f_64f(x);
		*yPointer++ = y;
	}
	return YepStatusOk;
}