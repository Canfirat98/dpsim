/* Copyright 2017-2021 Institute for Automation of Complex Power Systems,
 *                     EONERC, RWTH Aachen University
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at https://mozilla.org/MPL/2.0/.
 *********************************************************************************/

#include <dpsim-models/SP/SP_Ph1_SynchronGenerator5OrderVBR.h>

using namespace CPS;

SP::Ph1::SynchronGenerator5OrderVBR::SynchronGenerator5OrderVBR
    (const String & uid, const String & name, Logger::Level logLevel)
	: ReducedOrderSynchronGeneratorVBR(uid, name, logLevel),
	mEdq_t(mAttributes->create<Matrix>("Edq_t")),
	mEdq_s(mAttributes->create<Matrix>("Edq_s")) {

	mSGOrder = SGOrder::SG5Order;

	// model specific variables
	**mEdq_t = Matrix::Zero(2,1);
	**mEdq_s = Matrix::Zero(2,1);
	mEh_s = Matrix::Zero(2,1);
	mEh_t = Matrix::Zero(2,1);
}

SP::Ph1::SynchronGenerator5OrderVBR::SynchronGenerator5OrderVBR
	(const String & name, Logger::Level logLevel)
	: SynchronGenerator5OrderVBR(name, name, logLevel) {
}

void SP::Ph1::SynchronGenerator5OrderVBR::specificInitialization() {

	// initial voltage behind the transient reactance in the dq reference frame
	(**mEdq_t)(0,0) = 0.0;
	(**mEdq_t)(1,0) = (1 - mTaa / mTd0_t) * **mEf - (mLd - mLd_t - mYd) * (**mIdq)(0,0);

	// initial dq behind the subtransient reactance in the dq reference frame
	(**mEdq_s)(0,0) = (**mVdq)(0,0) - (mLq_s) * (**mIdq)(1,0);
	(**mEdq_s)(1,0) = (**mVdq)(1,0) + (mLd_s) * (**mIdq)(0,0);

	// initial history term behind the transient reactance
	mEh_t(0,0) = 0.0;
	mEh_t(1,0) = mAq_t * (**mIdq)(0,0) + mBq_t * (**mEdq_t)(1,0) + mDq_t * (**mEf) + mDq_t * mEf_prev;

	SPDLOG_LOGGER_DEBUG(mSLog,
		"\n--- Model specific initialization  ---"
		"\nSG model: 5th order"
		"\nInitial Ed_t (per unit): {:f}"
		"\nInitial Eq_t (per unit): {:f}"
		"\nInitial Ed_s (per unit): {:f}"
		"\nInitial Eq_s (per unit): {:f}"
		"\n--- Model specific initialization finished ---",
		(**mEdq_t)(0,0),
		(**mEdq_t)(1,0),
		(**mEdq_s)(0,0),
		(**mEdq_s)(1,0)
	);
	mSLog->flush();
}

void SP::Ph1::SynchronGenerator5OrderVBR::stepInPerUnit() {
	// calculate Edq_t at t=k
	(**mEdq_t)(0,0) = 0.0;
	(**mEdq_t)(1,0) = mAq_t * (**mIdq)(0,0) + mEh_t(1,0);

	// calculate Edq_s at t=k
	(**mEdq_s)(0,0) = (**mVdq)(0,0) - mLq_s * (**mIdq)(1,0);
	(**mEdq_s)(1,0) = (**mVdq)(1,0) + mLd_s * (**mIdq)(0,0);

	mDqToComplexA = get_DqToComplexATransformMatrix();
	mComplexAToDq = mDqToComplexA.transpose();

	// calculate resistance matrix at t=k+1
	calculateResistanceMatrix();

	// calculate history term behind the transient reactance
	mEh_t(0,0) = 0.0;
	mEh_t(1,0) = mAq_t * (**mIdq)(0,0) + mBq_t * (**mEdq_t)(1,0) + mDq_t * (**mEf) + mDq_t * mEf_prev;

	// calculate history term behind the subtransient reactance
	mEh_s(0,0) = mAd_s * (**mIdq)(1,0) + mCd_s * (**mEdq_s)(0,0);
	mEh_s(1,0) = mAq_s * (**mIdq)(0,0) + mBq_s * (**mEdq_t)(1,0) + mCq_s * (**mEdq_s)(1,0) + mDq_s * (**mEf) + mDq_s * mEf_prev;

	// convert Edq_t into the abc reference frame
	mEh_s = mDqToComplexA * mEh_s;
	mEvbr = Complex(mEh_s(0,0), mEh_s(1,0)) * mBase_V_RMS;
}
