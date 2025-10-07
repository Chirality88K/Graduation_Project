#ifndef FRENET_H
#define FRENET_H
#include "thirdparty/opennurbs/opennurbs.h"
#include "ChiralityLog.h"
#include <string>

class FrenetFrame
{
public:
	FrenetFrame()
	{
		mOrigin = ON_3dPoint(0, 0, 0);
		mAlpha = ON_3dVector(1, 0, 0);
		mBeta = ON_3dVector(0, 1, 0);
		mGamma = ON_3dVector(0, 0, 1);
	}
	FrenetFrame(ON_3dPoint o, ON_3dVector a, ON_3dVector b)
	{
		if (!Set(o, a, b))
		{
			CHIRALITY_ERROR(std::string("Frenet frame set incorrectly!!!!!!!"));
		}
	}
	bool Set(ON_3dPoint o, ON_3dVector a, ON_3dVector b)
	{
		mOrigin = o;
		a.Unitize();
		mAlpha = a;
		mBeta = b - ON_3dVector::DotProduct(a, b) * a;
		mBeta.Unitize();
		mGamma = ON_3dVector::CrossProduct(mAlpha, mBeta);
		return IsValid();
	}
	bool IsValid() const
	{
		return mAlpha.IsUnitVector() && mBeta.IsUnitVector() && mGamma.IsUnitVector() && mAlpha.IsPerpendicularTo(mBeta) &&
			   ((ON_3dVector::CrossProduct(mAlpha, mBeta) - mGamma).Length() < 1e-8);
	}
	ON_3dPoint GetPos()const { return mOrigin; }
	ON_3dVector GetAlpha() const { return mAlpha; }
	ON_3dVector GetBeta() const { return mBeta; }
	ON_3dVector GetGamma() const { return mGamma; }

private:
	ON_3dPoint mOrigin;
	ON_3dVector mAlpha;
	ON_3dVector mBeta;
	ON_3dVector mGamma;
};

#endif