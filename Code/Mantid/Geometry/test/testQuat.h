#ifndef MANTID_TESTQUAT__
#define MANTID_TESTQUAT__

#include <cxxtest/TestSuite.h>
#include <cmath>
#include <ostream>
#include "../inc/V3D.h"
#include "../inc/Quat.h" 

class testQuat : public CxxTest::TestSuite
{
private:

  Mantid::Geometry::Quat q,p;

public:
	void testoperatorbracket()
	{
		p[0]=0;
		p[1]=1;
		p[2]=2;
		p[3]=3;
		TS_ASSERT_EQUALS(p[0],0.0);
		TS_ASSERT_EQUALS(p[1],1.0);
		TS_ASSERT_EQUALS(p[2],2.0);
		TS_ASSERT_EQUALS(p[3],3.0);
	}
	void testEmptyConstructor()
	{
		TS_ASSERT_EQUALS(q[0],1.0);
		TS_ASSERT_EQUALS(q[1],0.0);
		TS_ASSERT_EQUALS(q[2],0.0);
		TS_ASSERT_EQUALS(q[3],0.0);
	}
	void testValueConstructor()
	{
		Mantid::Geometry::Quat q1(1,2,3,4);
		TS_ASSERT_EQUALS(q1[0],1.0);
		TS_ASSERT_EQUALS(q1[1],2.0);
		TS_ASSERT_EQUALS(q1[2],3.0);
		TS_ASSERT_EQUALS(q1[3],4.0);
	}
	void testAngleAxisConstructor()
	{
		Mantid::Geometry::V3D v(1,1,1);
		// Construct quaternion to represent rotation
		// of 45 degrees around the 111 axis.
		Mantid::Geometry::Quat q1(90.0,v);
		double c=1.0/sqrt(2.0);
		double s=c/sqrt(3.0);
		TS_ASSERT_DELTA(q1[0],c,0.000001);
		TS_ASSERT_DELTA(q1[1],s,0.000001);
		TS_ASSERT_DELTA(q1[2],s,0.000001);
		TS_ASSERT_DELTA(q1[3],s,0.000001);
	}
	void testoperatorassignmentfromdouble()
	{
		q(2,3,4,5);
		TS_ASSERT_EQUALS(q[0],2.0);
		TS_ASSERT_EQUALS(q[1],3.0);
		TS_ASSERT_EQUALS(q[2],4.0);
		TS_ASSERT_EQUALS(q[3],5.0);
	}
	void testoperatorassignmentfromangleaxis()
	{
		Mantid::Geometry::V3D v(1,1,1);
		q(90.0,v);
	        double c=1.0/sqrt(2.0);
		double s=c/sqrt(3.0);
		TS_ASSERT_DELTA(q[0],c,0.000001);
		TS_ASSERT_DELTA(q[1],s,0.000001);
		TS_ASSERT_DELTA(q[2],s,0.000001);
		TS_ASSERT_DELTA(q[3],s,0.000001);
	}
	void testoperatorequal()
	{
		q=p;
		TS_ASSERT_EQUALS(q[0],p[0]);
		TS_ASSERT_EQUALS(q[1],p[1]);
		TS_ASSERT_EQUALS(q[2],p[2]);
		TS_ASSERT_EQUALS(q[3],p[3]);
	}
	void testlenmethod()
	{
		q(1,2,3,4);
		TS_ASSERT_EQUALS(q.len(),sqrt(30.0));
	}
	void testlen2method()
	{
		q(1,2,3,4);
		TS_ASSERT_EQUALS(q.len2(),30.0);
	}
	void testinitmehtod()
	{
		q.init();
		TS_ASSERT_EQUALS(q[0],1);
		TS_ASSERT_EQUALS(q[1],0);
		TS_ASSERT_EQUALS(q[2],0);
		TS_ASSERT_EQUALS(q[3],0);
	}
	void testnormalizemethod()
	{
		q(2,2,2,2);
		q.normalize();
		TS_ASSERT_DELTA(q[0],0.125,0.000001);
		TS_ASSERT_DELTA(q[1],0.125,0.000001);
		TS_ASSERT_DELTA(q[2],0.125,0.000001)
		TS_ASSERT_DELTA(q[3],0.125,0.000001);
	}
	void testconjugatemethod()
	{
		q(1,1,1,1);
		q.conjugate();
		TS_ASSERT_EQUALS(q[0],1);
		TS_ASSERT_EQUALS(q[1],-1);
		TS_ASSERT_EQUALS(q[2],-1);
		TS_ASSERT_EQUALS(q[3],-1);
	}
	void testinversemethod()
	{
		q(2,3,4,5);
		Mantid::Geometry::Quat qinv(q);
		qinv.inverse();
		q*=qinv;
		TS_ASSERT_DELTA(q[0],1,0.000001);
		TS_ASSERT_DELTA(q[1],0,0.000001);
		TS_ASSERT_DELTA(q[2],0,0.000001);
		TS_ASSERT_DELTA(q[3],0,0.000001);
	}
	void testoperatorplus()
	{
		q(1,1,1,1);
		p(-1,2,1,3);
		Mantid::Geometry::Quat res;
		res=p+q;
		TS_ASSERT_EQUALS(res[0],0);
		TS_ASSERT_EQUALS(res[1],3);
		TS_ASSERT_EQUALS(res[2],2);
		TS_ASSERT_EQUALS(res[3],4);
	}
	void testoperatorminus()
	{
		q(1,1,1,1);
		p(-1,2,1,3);
		Mantid::Geometry::Quat res;
		res=p-q;
		TS_ASSERT_EQUALS(res[0],-2);
		TS_ASSERT_EQUALS(res[1],1);
		TS_ASSERT_EQUALS(res[2],0);
		TS_ASSERT_EQUALS(res[3],2);
	}
	void testoperatortimes()
	{
		q(1,1,1,1);
		p(-1,2,1,3);
		Mantid::Geometry::Quat res;
		res=p*q;
		TS_ASSERT_EQUALS(res[0],-7);
		TS_ASSERT_EQUALS(res[1],-1);
		TS_ASSERT_EQUALS(res[2],1);
		TS_ASSERT_EQUALS(res[3],3);
	}
	void testoperatordoublequal()
	{
		p=q;
		TS_ASSERT(p==q);
		q(1,4,5,6);
		TS_ASSERT(p!=q);
	}
	void testoperatornotequal()
	{
		q(1,2,3,4);
		TS_ASSERT(p!=q);
		p=q;
		TS_ASSERT(!(p!=q));
	}
	};

#endif