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


extern "C" YEP_LOCAL_SYMBOL YepStatus _yepMath_EvaluatePolynomial_V32fV32f_V32f_Default(const Yep32f *YEP_RESTRICT coefPointer, const Yep32f *YEP_RESTRICT xPointer, Yep32f *YEP_RESTRICT yPointer, YepSize coefCount, YepSize length) {
	if YEP_UNLIKELY(coefPointer == YEP_NULL_POINTER) {
		return YepStatusNullPointer;
	}
	if YEP_UNLIKELY(yepBuiltin_GetPointerMisalignment(coefPointer, sizeof(Yep32f)) != 0) {
		return YepStatusMisalignedPointer;
	}
	if YEP_UNLIKELY(xPointer == YEP_NULL_POINTER) {
		return YepStatusNullPointer;
	}
	if YEP_UNLIKELY(yepBuiltin_GetPointerMisalignment(xPointer, sizeof(Yep32f)) != 0) {
		return YepStatusMisalignedPointer;
	}
	if YEP_UNLIKELY(yPointer == YEP_NULL_POINTER) {
		return YepStatusNullPointer;
	}
	if YEP_UNLIKELY(yepBuiltin_GetPointerMisalignment(yPointer, sizeof(Yep32f)) != 0) {
		return YepStatusMisalignedPointer;
	}
	if YEP_UNLIKELY(coefCount == 0) {
		return YepStatusInvalidArgument;
	}
	if YEP_UNLIKELY(coefCount == 0) {
		return YepStatusInvalidArgument;
	}
	while (length-- != 0) {
		const Yep32f x = *xPointer++;
		Yep32f y = coefPointer[coefCount - 1];
		for (YepSize coefIndex = coefCount - 1; coefIndex != 0; coefIndex--) {
			const Yep32f coef = coefPointer[coefIndex - 1];
			y = yepBuiltin_MultiplyAdd_32f32f32f_32f(y, x, coef);
		}
		*yPointer++ = y;
	}
	return YepStatusOk;
}

extern "C" YEP_LOCAL_SYMBOL YepStatus _yepMath_EvaluatePolynomial_V64fV64f_V64f_Default(const Yep64f *YEP_RESTRICT coefPointer, const Yep64f *YEP_RESTRICT xPointer, Yep64f *YEP_RESTRICT yPointer, YepSize coefCount, YepSize length) {
	if YEP_UNLIKELY(coefPointer == YEP_NULL_POINTER) {
		return YepStatusNullPointer;
	}
	if YEP_UNLIKELY(yepBuiltin_GetPointerMisalignment(coefPointer, sizeof(Yep64f)) != 0) {
		return YepStatusMisalignedPointer;
	}
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
	if YEP_UNLIKELY(coefCount == 0) {
		return YepStatusInvalidArgument;
	}
	if YEP_UNLIKELY(coefCount == 0) {
		return YepStatusInvalidArgument;
	}
	while (length-- != 0) {
		const Yep64f x = *xPointer++;
		Yep64f y = coefPointer[coefCount - 1];
		for (YepSize coefIndex = coefCount - 1; coefIndex != 0; coefIndex--) {
			const Yep64f coef = coefPointer[coefIndex - 1];
			y = yepBuiltin_MultiplyAdd_64f64f64f_64f(y, x, coef);
		}
		*yPointer++ = y;
	}
	return YepStatusOk;
}
