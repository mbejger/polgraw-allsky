/*
 *                            Yeppp! library header
 *                   This file is auto-generated by Peach-Py,
 *        Portable Efficient Assembly Code-generator in Higher-level Python,
 *                    part of the Yeppp! library infrastrure
 * 
 * This file is part of Yeppp! library and licensed under the New BSD license.
 * 
 * Copyright (C) 2010-2012 Marat Dukhan
 * Copyright (C) 2012-2013 Georgia Institute of Technology
 * All rights reserved.
 *  
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *     * Redistributions of source code must retain the above copyright
 *       notice, this list of conditions and the following disclaimer.
 *     * Redistributions in binary form must reproduce the above copyright
 *       notice, this list of conditions and the following disclaimer in the
 *       documentation and/or other materials provided with the distribution.
 *     * Neither the name of the Georgia Institute of Technology nor the
 *       names of its contributors may be used to endorse or promote products
 *       derived from this software without specific prior written permission.
 * 
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
 * ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
 * WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
 * DISCLAIMED. IN NO EVENT SHALL <COPYRIGHT HOLDER> BE LIABLE FOR ANY
 * DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
 * (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
 * LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
 * ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 */

#pragma once

#include <yepPredefines.h>
#include <yepTypes.h>

#ifdef __cplusplus
	extern "C" {
#endif

/** @defgroup yepMath yepMath.h: vector mathematical functions. */

/**
 * @ingroup yepMath
 * @defgroup yepMath_Log	Natural logarithm
 */

/**
 * @ingroup	yepMath_Log
 * @brief	Computes natural logarithm on an array of double precision (64-bit) floating-point elements.
 * @param[in]	x	Pointer to the array of elements on which logarithm will be computed.
 * @param[out]	y	Pointer the array where the computed logarithms will be stored.
 * @param[in]	length	Length of the arrays specified by @a x and @a y.
 * @retval	#YepStatusOk	The computation finished successfully.
 * @retval	#YepStatusNullPointer	@a x or @a y argument is null.
 * @retval	#YepStatusMisalignedPointer	@a x or @a y argument is not naturally aligned.
 * @par	Optimized implementations
 *     		<table>
 *     			<tr><th>Architecture</th><th>Target microarchitecture</th><th>Required instruction extensions</th></tr>
 *     			<tr><td>x86-64</td><td>Intel Nehalem</td><td>SSE, SSE2, SSE4.1</td></tr>
 *     			<tr><td>x86-64</td><td>Intel Sandy Bridge</td><td>AVX</td></tr>
 *     			<tr><td>x86-64</td><td>Intel Haswell</td><td>AVX, AVX2, FMA3</td></tr>
 *     			<tr><td>x86-64</td><td>AMD K10</td><td>SSE, SSE2</td></tr>
 *     			<tr><td>x86-64</td><td>AMD Bulldozer</td><td>AVX, FMA4, XOP</td></tr>
 *     			<tr><td>x86-64</td><td>AMD Bobcat</td><td>SSE, SSE2</td></tr>
 *     		</table>
 */
YEP_PUBLIC_SYMBOL enum YepStatus YEPABI yepMath_Log_V64f_V64f(const Yep64f *YEP_RESTRICT x, Yep64f *YEP_RESTRICT y, YepSize length);
/**
 * @ingroup yepMath
 * @defgroup yepMath_Exp	Base-e exponent
 */

/**
 * @ingroup	yepMath_Exp
 * @brief	Computes base-e exponent on an array of double precision (64-bit) floating-point elements.
 * @param[in]	x	Pointer to the array of elements on which exponent will be computed.
 * @param[out]	y	Pointer the array where the computed exponents will be stored.
 * @param[in]	length	Length of the arrays specified by @a x and @a y.
 * @retval	#YepStatusOk	The computation finished successfully.
 * @retval	#YepStatusNullPointer	@a x or @a y argument is null.
 * @retval	#YepStatusMisalignedPointer	@a x or @a y argument is not naturally aligned.
 * @par	Optimized implementations
 *     		<table>
 *     			<tr><th>Architecture</th><th>Target microarchitecture</th><th>Required instruction extensions</th></tr>
 *     			<tr><td>x86-64</td><td>Intel Nehalem</td><td>SSE, SSE2, SSE4.1</td></tr>
 *     			<tr><td>x86-64</td><td>Intel Sandy Bridge</td><td>AVX</td></tr>
 *     			<tr><td>x86-64</td><td>Intel Haswell</td><td>AVX, AVX2, FMA3</td></tr>
 *     			<tr><td>x86-64</td><td>AMD K10</td><td>CMOV, SSE, SSE2</td></tr>
 *     			<tr><td>x86-64</td><td>AMD Bulldozer</td><td>AVX, FMA4</td></tr>
 *     			<tr><td>x86-64</td><td>AMD Bobcat</td><td>CMOV, SSE, SSE2</td></tr>
 *     		</table>
 */
YEP_PUBLIC_SYMBOL enum YepStatus YEPABI yepMath_Exp_V64f_V64f(const Yep64f *YEP_RESTRICT x, Yep64f *YEP_RESTRICT y, YepSize length);
/**
 * @ingroup yepMath
 * @defgroup yepMath_Sin	Sine
 */

/**
 * @ingroup	yepMath_Sin
 * @brief	Computes sine on an array of double precision (64-bit) floating-point elements.
 * @param[in]	x	Pointer to the array of elements on which sine will be computed.
 * @param[out]	y	Pointer the array where the computed sines will be stored.
 * @param[in]	length	Length of the arrays specified by @a x and @a y.
 * @retval	#YepStatusOk	The computation finished successfully.
 * @retval	#YepStatusNullPointer	@a x or @a y argument is null.
 * @retval	#YepStatusMisalignedPointer	@a x or @a y argument is not naturally aligned.
 * @par	Optimized implementations
 *     		<table>
 *     			<tr><th>Architecture</th><th>Target microarchitecture</th><th>Required instruction extensions</th></tr>
 *     			<tr><td>x86-64</td><td>Intel Nehalem</td><td>SSE, SSE2, SSE4.1</td></tr>
 *     			<tr><td>x86-64</td><td>Intel Sandy Bridge</td><td>AVX</td></tr>
 *     			<tr><td>x86-64</td><td>Intel Haswell</td><td>AVX, AVX2, FMA3</td></tr>
 *     			<tr><td>x86-64</td><td>AMD Bulldozer</td><td>AVX, FMA4</td></tr>
 *     		</table>
 */
YEP_PUBLIC_SYMBOL enum YepStatus YEPABI yepMath_Sin_V64f_V64f(const Yep64f *YEP_RESTRICT x, Yep64f *YEP_RESTRICT y, YepSize length);
/**
 * @ingroup yepMath
 * @defgroup yepMath_Cos	Cosine
 */

/**
 * @ingroup	yepMath_Cos
 * @brief	Computes cosine on an array of double precision (64-bit) floating-point elements.
 * @param[in]	x	Pointer to the array of elements on which cosine will be computed.
 * @param[out]	y	Pointer the array where the computed cosines will be stored.
 * @param[in]	length	Length of the arrays specified by @a x and @a y.
 * @retval	#YepStatusOk	The computation finished successfully.
 * @retval	#YepStatusNullPointer	@a x or @a y argument is null.
 * @retval	#YepStatusMisalignedPointer	@a x or @a y argument is not naturally aligned.
 * @par	Optimized implementations
 *     		<table>
 *     			<tr><th>Architecture</th><th>Target microarchitecture</th><th>Required instruction extensions</th></tr>
 *     			<tr><td>x86-64</td><td>Intel Nehalem</td><td>SSE, SSE2, SSE4.1</td></tr>
 *     			<tr><td>x86-64</td><td>Intel Sandy Bridge</td><td>AVX</td></tr>
 *     			<tr><td>x86-64</td><td>Intel Haswell</td><td>AVX, AVX2, FMA3</td></tr>
 *     			<tr><td>x86-64</td><td>AMD Bulldozer</td><td>AVX, FMA4</td></tr>
 *     		</table>
 */
YEP_PUBLIC_SYMBOL enum YepStatus YEPABI yepMath_Cos_V64f_V64f(const Yep64f *YEP_RESTRICT x, Yep64f *YEP_RESTRICT y, YepSize length);
/**
 * @ingroup yepMath
 * @defgroup yepMath_Tan	Tangent
 */

/**
 * @ingroup	yepMath_Tan
 * @brief	Computes tangent on an array of double precision (64-bit) floating-point elements.
 * @param[in]	x	Pointer to the array of elements on which tangent will be computed.
 * @param[out]	y	Pointer the array where the computed tangents will be stored.
 * @param[in]	length	Length of the arrays specified by @a x and @a y.
 * @retval	#YepStatusOk	The computation finished successfully.
 * @retval	#YepStatusNullPointer	@a x or @a y argument is null.
 * @retval	#YepStatusMisalignedPointer	@a x or @a y argument is not naturally aligned.
 * @par	Optimized implementations
 *     		<table>
 *     			<tr><th>Architecture</th><th>Target microarchitecture</th><th>Required instruction extensions</th></tr>
 *     			<tr><td>x86-64</td><td>AMD Bulldozer</td><td>AVX, FMA4</td></tr>
 *     		</table>
 */
YEP_PUBLIC_SYMBOL enum YepStatus YEPABI yepMath_Tan_V64f_V64f(const Yep64f *YEP_RESTRICT x, Yep64f *YEP_RESTRICT y, YepSize length);
/**
 * @ingroup yepMath
 * @defgroup yepMath_EvaluatePolynomial	Polynomial evaluation
 */

/**
 * @ingroup	yepMath_EvaluatePolynomial
 * @brief	Evaluates polynomial with single precision (32-bit) floating-point coefficients on an array of single precision (32-bit) floating-point elements.
 * @param[in]	x	Pointer to the array of elements on which the polynomial will be evaluated.
 * @param[in]	coef	Pointer to the array of polynomial coefficients.
 * @param[out]	y	Pointer the array where the result of polynomial evaluation will be stored.
 * @param[in]	coefCount	Number of polynomial coefficients. Should equal the polynomial degree plus one.
 * @param[in]	length	Length of the arrays specified by @a x and @a y.
 * @retval	#YepStatusOk	The computation finished successfully.
 * @retval	#YepStatusNullPointer	@a coef, @a x or @a y argument is null.
 * @retval	#YepStatusMisalignedPointer	@a coef, @a x or @a y argument is not naturally aligned.
 * @retval	#YepStatusInvalidArgument	@a coefCount argument is zero.
 * @par	Optimized implementations
 *     		<table>
 *     			<tr><th>Architecture</th><th>Target microarchitecture</th><th>Required instruction extensions</th></tr>
 *     			<tr><td>x86-64</td><td>Intel Nehalem</td><td>SSE</td></tr>
 *     			<tr><td>x86-64</td><td>Intel Sandy Bridge</td><td>AVX</td></tr>
 *     			<tr><td>x86-64</td><td>Intel Haswell</td><td>SSE, AVX, FMA3</td></tr>
 *     			<tr><td>x86-64</td><td>Intel Bonnell</td><td>SSE</td></tr>
 *     			<tr><td>x86-64</td><td>AMD Bulldozer</td><td>AVX, FMA4</td></tr>
 *     			<tr><td>ARM</td><td>ARM Cortex-A9</td><td>VFP2, NEON</td></tr>
 *     		</table>
 */
YEP_PUBLIC_SYMBOL enum YepStatus YEPABI yepMath_EvaluatePolynomial_V32fV32f_V32f(const Yep32f *YEP_RESTRICT coef, const Yep32f *YEP_RESTRICT x, Yep32f *YEP_RESTRICT y, YepSize coefCount, YepSize length);
/**
 * @ingroup	yepMath_EvaluatePolynomial
 * @brief	Evaluates polynomial with double precision (64-bit) floating-point coefficients on an array of double precision (64-bit) floating-point elements.
 * @param[in]	x	Pointer to the array of elements on which the polynomial will be evaluated.
 * @param[in]	coef	Pointer to the array of polynomial coefficients.
 * @param[out]	y	Pointer the array where the result of polynomial evaluation will be stored.
 * @param[in]	coefCount	Number of polynomial coefficients. Should equal the polynomial degree plus one.
 * @param[in]	length	Length of the arrays specified by @a x and @a y.
 * @retval	#YepStatusOk	The computation finished successfully.
 * @retval	#YepStatusNullPointer	@a coef, @a x or @a y argument is null.
 * @retval	#YepStatusMisalignedPointer	@a coef, @a x or @a y argument is not naturally aligned.
 * @retval	#YepStatusInvalidArgument	@a coefCount argument is zero.
 * @par	Optimized implementations
 *     		<table>
 *     			<tr><th>Architecture</th><th>Target microarchitecture</th><th>Required instruction extensions</th></tr>
 *     			<tr><td>x86-64</td><td>Intel Nehalem</td><td>SSE2, SSE3</td></tr>
 *     			<tr><td>x86-64</td><td>Intel Sandy Bridge</td><td>AVX</td></tr>
 *     			<tr><td>x86-64</td><td>Intel Haswell</td><td>SSE, AVX, FMA3</td></tr>
 *     			<tr><td>x86-64</td><td>Intel Bonnell</td><td>SSE, SSE2</td></tr>
 *     			<tr><td>x86-64</td><td>AMD Bulldozer</td><td>AVX, FMA4</td></tr>
 *     			<tr><td>ARM</td><td>ARM Cortex-A9</td><td>VFP2, VFPd32</td></tr>
 *     		</table>
 */
YEP_PUBLIC_SYMBOL enum YepStatus YEPABI yepMath_EvaluatePolynomial_V64fV64f_V64f(const Yep64f *YEP_RESTRICT coef, const Yep64f *YEP_RESTRICT x, Yep64f *YEP_RESTRICT y, YepSize coefCount, YepSize length);
#ifdef __cplusplus
	} // extern "C"
#endif