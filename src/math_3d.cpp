/*
(usage)
math.exe mode arg1 arg2 ...

mode       必須
arg1...N   modeに応じて必要な個数だけ


(行列の並び)
右手系,row vectors.
*3x3
	m00, m01, m02 | Tx
	m10, m11, m12 | Ty
	m20, m21, m22 | Tz
	--------------+----
	              | 1.0

(参考資料)
Graphics gems
Game programming gems
Game developer magazine
Java 3d
Real time rendering
Pysics rendering
ゲームプログラマーになる前に呼んでおくこと
クォータニオン入門　金谷一郎
CEDEC講演資料
http://www.euclideanspace.com/maths/geometry/rotations/conversions/matrixToAngle/index.htm
グーグル先生

（角度の単位）
Dcc-Tool(3dsMAX/MAYA/XSI)の画面表示とデザインとの意思疎通のため
ラジアンより度数表記を優先しています。
*/

extern "C"{
	#include "gems_iv/EulerAngles.h"	/* graphics gems iv*/
};

#include "stdafx.h"

#define EPSILON		0.00001f
#define	PI			3.14159265358979323
typedef unsigned int 		uint32;
typedef unsigned __int64 	uint64;

static void usage(){
    printf(
        "math.exe mode arg1 arg2 ...\n"
        "mode\n"
        "  f32_hex      float             -> hex(uint32)\n"
        "  f64_hex      double            -> hex(uint64)\n"
        "  hex_f32      hex(uint32)       -> float\n"
        "  hex_f64      hex(uint64)       -> double\n"
        "  2_10         bin               -> decimal\n"
        "  2_16         bin               -> hex\n"
        "  10_2         decimal           -> bin\n"
        "  10_16        decimal           -> hex\n"
        "  16_2         hex               -> bin\n"
        "  16_10        hex               -> decimal\n"
        "  r_d          radian            -> degree\n"
        "  d_r          degree            -> radian\n"
        "  v_h          vertical fov(degree)    -> horizontal fov(degree)\n"
        "  h_v          horizontal fov(degree)  -> vertical fov(degree)\n"
        "  ql           quaternion length\n"
        "  qn           quaternion normalize\n"
        "  qi           quaternion invert\n"
        "  q_m33        quaternion        -> matrix33\n"
        "  q_a          quaternion        -> axis angle\n"
        "  q_xyz        quaternion        -> euler xyz(degree)\n"
        "  x_m33        x(degree)         -> matrix33\n"
        "  y_m33        y(degree)         -> matrix33\n"
        "  z_m33        z(degree)         -> matrix33\n"
        "  yxz_m33      euler yxz(degree) -> matrix33\n"
        "  zxy_m33      euler zxy(degree) -> matrix33\n"
        "  zyx_m33      euler zyx(degree) -> matrix33\n"
        "  yzx_m33      euler yzx(degree) -> matrix33\n"
        "  xzy_m33      euler xzy(degree) -> matrix33\n"
        "  xyz_m33      euler xyz(degree) -> matrix33\n"
        "  m33t         matrix33 transpose\n"
        "  m33d         matrix33 determinant\n"
        "  m33i         matrix33 invert\n"
        "  m33_q        matrix33          -> quaternion\n"
        "  m33_a        matrix33          -> axis angle\n"
        "  m33_s        matrix33          -> scale xyz\n"
        "  m33_yxz      matrix33          -> yxz(degree)\n"
        "  m33_zxy      matrix33          -> zxy(degree)\n"
        "  m33_zyx      matrix33          -> zyx(degree)\n"
        "  m33_yzx      matrix33          -> yzx(degree)\n"
        "  m33_xzy      matrix33          -> xzy(degree)\n"
        "  m33_xyz      matrix33          -> xyz(degree)\n"
        "\n"
        "unit\n"
        "  degree      0-360\n"
        "  radian      0-2pi\n"
        "  quaternion  x y z w\n"
        "  axis_angle  x y z w (xyz=axis w=angle)\n"
        "  matrix33    m00 m01 m02|Tx\n"
        "              m10 m11 m12|Ty\n"
        "              m20 m21 m22|Tz\n"
        "              -----------+----\n"
        "              0.0 0.0 0.0|1.0\n"
        "ex.\n"
        "math f32_hex 1.0\n"
        " -> 3f800000\n"
        "math hex_f32 3f800000\n"
        " -> 1.000000\n"
        "math d_r 180\n"
        " -> 3.141593\n"
        "math 16_10 a 0xb ffh\n"
        " -> 10 11 255\n"
        "math 10_16 64 255 -1\n"
        " ->0x40 0xff 0xffffffffffffffff\n"
        "math m33t 1 2 3 4 5 6 7 8 9\n"
        " -> 1.000000 4.000000 7.000000\n"
        "    2.000000 5.000000 8.000000\n"
        "    3.000000 6.000000 9.000000\n"
    );
}

typedef std::deque<std::string>		TString1D;
struct TArgInfo{
	std::string		m_mode;
	TString1D		m_arg;
};


static float deg2rad(float d){
	return d * (PI / 180.f);
}
static float rad2deg(float r){
	return r * (180.f / PI);
}
static double deg2rad(double d){
	return d * (PI / 180.0);
}
static double rad2deg(double r){
	return r * (180.0 / PI);
}

//回転行列(3x3)
struct TMatrix{
	TMatrix(){
		m00=1.f; m01=0.f; m02=0.f;
		m10=0.f; m11=1.f; m12=0.f;
		m20=0.f; m21=0.f; m22=1.f;
	};
	TMatrix(const TString1D &in){
		m00=atof(in[0].c_str()); m01=atof(in[1].c_str()); m02=atof(in[2].c_str());
		m10=atof(in[3].c_str()); m11=atof(in[4].c_str()); m12=atof(in[5].c_str());
		m20=atof(in[6].c_str()); m21=atof(in[7].c_str()); m22=atof(in[8].c_str());

	};

	TMatrix(float n00, float n01, float n02,
            float n10, float n11, float n12,
            float n20, float n21, float n22)
	:	m00(n00), m01(n01), m02(n02),
    	m10(n10), m11(n11), m12(n12),
    	m20(n20), m21(n21), m22(n22)
    {
    }
	void Set(float m00, float m01, float m02,
             float m10, float m11, float m12,
             float m20, float m21, float m22) {
    	this->m00 = m00; this->m01 = m01; this->m02 = m02;
    	this->m10 = m10; this->m11 = m11; this->m12 = m12;
    	this->m20 = m20; this->m21 = m21; this->m22 = m22;
	}


    void Invert(){
    	float s = 1.f/Determinant();	//oops! div 0;
    	Set(
        	m11*m22 - m12*m21, m02*m21 - m01*m22, m01*m12 - m02*m11,
        	m12*m20 - m10*m22, m00*m22 - m02*m20, m02*m10 - m00*m12,
        	m10*m21 - m11*m20, m01*m20 - m00*m21, m00*m11 - m01*m10
        );
    	Mul(s);
	};

	float Determinant()const{
		return 	  m00*(m11*m22 - m21*m12)
        		- m01*(m10*m22 - m20*m12)
        		+ m02*(m10*m21 - m20*m11);
	}

	void Mul(const TMatrix& m1, const TMatrix& m2) {
	    Set(
	        m1.m00*m2.m00 + m1.m01*m2.m10 + m1.m02*m2.m20,
	        m1.m00*m2.m01 + m1.m01*m2.m11 + m1.m02*m2.m21,
	        m1.m00*m2.m02 + m1.m01*m2.m12 + m1.m02*m2.m22,

	        m1.m10*m2.m00 + m1.m11*m2.m10 + m1.m12*m2.m20,
	        m1.m10*m2.m01 + m1.m11*m2.m11 + m1.m12*m2.m21,
	        m1.m10*m2.m02 + m1.m11*m2.m12 + m1.m12*m2.m22,

	        m1.m20*m2.m00 + m1.m21*m2.m10 + m1.m22*m2.m20,
	        m1.m20*m2.m01 + m1.m21*m2.m11 + m1.m22*m2.m21,
	        m1.m20*m2.m02 + m1.m21*m2.m12 + m1.m22*m2.m22
	    );
	}

	void Mul(const TMatrix& m1) {
    	Mul(*this, m1);
	}

	void Mul(float scalar) {
    	m00 *= scalar; m01 *= scalar;  m02 *= scalar;
    	m10 *= scalar; m11 *= scalar;  m12 *= scalar;
    	m20 *= scalar; m21 *= scalar;  m22 *= scalar;
    }

    void Transpose() {
	    float tmp = m01;
	    m01 = m10;
	    m10 = tmp;

	    tmp = m02;
	    m02 = m20;
	    m20 = tmp;

	    tmp = m12;
	    m12 = m21;
	    m21 = tmp;
	}

	void RotX(float angle) {
		float c = cos(angle);
		float s = sin(angle);
	    m00 = 1.0; m01 = 0.0; m02 = 0.0;
	    m10 = 0.0; m11 = c;   m12 = -s;
	    m20 = 0.0; m21 = s;   m22 = c;
	}

	void RotY(float angle) {
		float c = cos(angle);
		float s = sin(angle);
	    m00 = c;   m01 = 0.0; m02 = s;
	    m10 = 0.0; m11 = 1.0; m12 = 0.0;
	    m20 = -s;  m21 = 0.0; m22 = c;
	}

	void RotZ(float angle) {
		float c = cos(angle);
		float s = sin(angle);
	    m00 = c;   m01 = -s;  m02 = 0.0;
	    m10 = s;   m11 = c;   m12 = 0.0;
	    m20 = 0.0; m21 = 0.0; m22 = 1.0;
	}

	void ExtractScale(float &x, float &y, float &z)const{
		x = sqrt(m00*m00 + m01*m01 + m02*m02);
		y = sqrt(m10*m10 + m11*m11 + m12*m12);
		z = sqrt(m20*m20 + m21*m21 + m22*m22);
	}

	TMatrix& operator*=(const TMatrix& m1) {
        Mul(m1);
        return *this;
    }

	union{
		struct {
			float m00, m01, m02;
			float m10, m11, m12;
			float m20, m21, m22;
		};
		float m[3][3];
	};
};

TMatrix operator*(const TMatrix& m1, const TMatrix& m2) {
    return (TMatrix(m1)).operator*=(m2);
}

struct CQuaternion{
	CQuaternion(){
		x=0.f; y=0.f; z=0.f; w=0.f;
	};
	CQuaternion(float _x, float _y, float _z, float _w){
		x=_x;
		y=_y;
		z=_z;
		w=_w;
	};
	CQuaternion(const std::string &_x, const std::string &_y, const std::string &_z, const std::string &_w){
		x = atof(_x.c_str());
		y = atof(_y.c_str());
		z = atof(_z.c_str());
		w = atof(_w.c_str());
	}

	float Square()const{
		return x*x + y*y + z*z + w*w;
	};
	float Length()const{
		return sqrtf(Square());
	};

	float x,y,z,w;

};

struct TTuple3{
	TTuple3(){
		x=0.f;
		y=0.f;
		z=0.f;
	};


	float x,y,z;
};

//memo:微妙なtypedefですが・・・
typedef CQuaternion	TAxisAngle;


template<typename FUNC>
static bool MatrixSub(const TArgInfo &info, FUNC f, size_t unit){
	if(! (unit<=info.m_arg.size())){
		return false;
	}

	TArgInfo 	tmp = info;
	while(unit <= tmp.m_arg.size()){
		f(tmp);
		tmp.m_arg.erase(tmp.m_arg.begin(), tmp.m_arg.begin()+unit);
	}
	if(tmp.m_arg.size()){
		printf("Warning:Remain %d.\n", tmp.m_arg.size());
	}
	return true;
}

static void ConvertMatrix(HMatrix&out, const TMatrix&in){
	TMatrix tmp = in;
	//tmp.Transpose();

	out[0][0]=tmp.m[0][0];
	out[0][1]=tmp.m[0][1];
	out[0][2]=tmp.m[0][2];
	out[0][3]=0.f;

	out[1][0]=tmp.m[1][0];
	out[1][1]=tmp.m[1][1];
	out[1][2]=tmp.m[1][2];
	out[1][3]=0.f;

	out[2][0]=tmp.m[2][0];
	out[2][1]=tmp.m[2][1];
	out[2][2]=tmp.m[2][2];
	out[2][3]=0.f;

	out[3][0]=0.f;
	out[3][1]=0.f;
	out[3][2]=0.f;
	out[3][3]=1.f;
}

static void ConvertMatrix(TMatrix&out, const HMatrix&in){
	out.m[0][0]=in[0][0];
	out.m[0][1]=in[0][1];
	out.m[0][2]=in[0][2];

	out.m[1][0]=in[1][0];
	out.m[1][1]=in[1][1];
	out.m[1][2]=in[1][2];

	out.m[2][0]=in[2][0];
	out.m[2][1]=in[2][1];
	out.m[2][2]=in[2][2];

	//out.Transpose();
}

static void QuatToMatMain( TMatrix &out, const CQuaternion &in){
	float x = in.x;
	float y = in.y;
	float z = in.z;
	float w = in.w;

    float x2 = x*2.f,  y2 = y*2.f,  z2 = z*2.f;
    float wx = w*x2, wy = w*y2, wz = w*z2;
    float xx = x*x2, xy = x*y2, xz = x*z2;
    float yy = y*y2, yz = y*z2, zz = z*z2;

    out.m00 = 1.0 - (yy + zz); out.m01 = xy - wz;         out.m02 = xz + wy;
    out.m10 = xy + wz;         out.m11 = 1.0 - (xx + zz); out.m12 = yz - wx;
    out.m20 = xz - wy;         out.m21 = yz + wx;         out.m22 = 1.0 - (xx + yy);
}

static void QuatToAxisAngleMain(TAxisAngle&out, const CQuaternion&in){
	float angle = 2.0f * acos(in.w);
	float ax	= sin( 0.5f * angle );

	if(fabs(ax) < EPSILON){
		out.x = 1.f;
		out.y = 0.f;
		out.z = 0.f;
	}else{
		float scl = 1.f/ax;
		out.x = in.x * scl;
		out.y = in.y * scl;
		out.z = in.z * scl;
	}
	out.w = angle;
}

static void Mat33ToQuatMain(CQuaternion &out, const TMatrix &in){
	float s=0.f;
	float tr = in.m00 + in.m11 + in.m22;
	if (tr >= 0.0) {
	    s = sqrt(tr + 1.0);
	    out.w = s*0.5;
	    s = 0.5/s;
	    out.x = (in.m21 - in.m12)*s;
	    out.y = (in.m02 - in.m20)*s;
	    out.z = (in.m10 - in.m01)*s;
	} else {
	    float maxm = std::max(in.m00, std::max(in.m11, in.m22));
	    if (maxm == in.m00) {
            s = sqrt(in.m00 - (in.m11 + in.m22) + 1.0);
            out.x = s*0.5;
            s = 0.5/s;
            out.y = (in.m01 + in.m10)*s;
            out.z = (in.m20 + in.m02)*s;
            out.w = (in.m21 - in.m12)*s;
	    } else if (maxm == in.m11) {
            s = sqrt(in.m11 - (in.m22 + in.m00) + 1.0);
            out.y = s*0.5;
            s = 0.5/s;
            out.z = (in.m12 + in.m21)*s;
            out.x = (in.m01 + in.m10)*s;
            out.w = (in.m02 - in.m20)*s;
	    } else {
            s = sqrt(in.m22 - (in.m00 + in.m11) + 1.0);
            out.z = s*0.5;
            s = 0.5/s;
            out.x = (in.m20 + in.m02)*s;
            out.y = (in.m12 + in.m21)*s;
            out.w = (in.m10 - in.m01)*s;
	    }
	}
}

static void Mat33ToAAngMain(TAxisAngle&out, const TMatrix &in){
#if 1
	double x=0.0;
	double y=0.0;
	double z=0.0;
	double angle=0.0;
	double epsilon = 0.01; // margin to allow for rounding errors
	double epsilon2 = 0.1; // margin to distinguish between 0 and 180 degrees

	if ( (fabs(in.m01-in.m10)< epsilon)
	  && (fabs(in.m02-in.m20)< epsilon)
	  && (fabs(in.m12-in.m21)< epsilon)) {
		// singularity found
		// first check for identity matrix which must have +1 for all terms
		//  in leading diagonaland zero in other terms
		if ( (fabs(in.m01+in.m10) < epsilon2)
		  && (fabs(in.m02+in.m20) < epsilon2)
		  && (fabs(in.m12+in.m21) < epsilon2)
		  && (fabs(in.m00+in.m11+in.m22-3) < epsilon2)) {
			// this singularity is identity matrix so angle = 0
			out.x = 1.0;
			out.y = 0.0;
			out.z = 0.0;
			out.w = 0.0;
			return ;
			//return new axisAngle(0,1,0,0); // zero angle, arbitrary axis
		}
		// otherwise this singularity is angle = 180
		angle = PI;
		double xx = (in.m00+1)/2;
		double yy = (in.m11+1)/2;
		double zz = (in.m22+1)/2;
		double xy = (in.m01+in.m10)/4;
		double xz = (in.m02+in.m20)/4;
		double yz = (in.m12+in.m21)/4;
		if ((xx > yy) && (xx > zz)) { // m[0][0] is the largest diagonal term
			if (xx< epsilon) {
				x = 0;
				y = 0.7071;
				z = 0.7071;
			} else {
				x = sqrt(xx);
				y = xy/x;
				z = xz/x;
			}
		} else if (yy > zz) { // m[1][1] is the largest diagonal term
			if (yy< epsilon) {
				x = 0.7071;
				y = 0;
				z = 0.7071;
			} else {
				y = sqrt(yy);
				x = xy/y;
				z = yz/y;
			}
		} else { // m22 is the largest diagonal term so base result on this
			if (zz< epsilon) {
				x = 0.7071;
				y = 0.7071;
				z = 0;
			} else {
				z = sqrt(zz);
				x = xz/z;
				y = yz/z;
			}
		}
		out.x = x;
		out.y = y;
		out.z = z;
		out.w = angle;
		//return new axisAngle(angle,x,y,z); // return 180 deg rotation
	}

	// as we have reached here there are no singularities so we can handle normally
	double s = sqrt(( in.m21 - in.m12)*(in.m21 - in.m12)
					+(in.m02 - in.m20)*(in.m02 - in.m20)
					+(in.m10 - in.m01)*(in.m10 - in.m01)); // used to normalise
	if (fabs(s) < 0.001) s=1;
		// prevent divide by zero, should not happen if matrix is orthogonal and should be
		// caught by singularity test above, but I've left it in just in case
	angle = acos(( in.m00 + in.m11 + in.m22 - 1)/2);
	x = (in.m21 - in.m12)/s;
	y = (in.m02 - in.m20)/s;
	z = (in.m10 - in.m01)/s;
	out.x = x;
	out.y = y;
	out.z = z;
	out.w = angle;
   //return new axisAngle(angle,x,y,z);
#else
	float cos_val = (in.m00 + in.m11 + in.m22 - 1.0)*0.5;
	float x = in.m21 - in.m12;
	float y = in.m02 - in.m20;
	float z = in.m10 - in.m01;
	float sin_val = 0.5*sqrt(x*x + y*y + z*z);
	float w = atan2(sin_val, cos_val);

	out.x = x;
	out.y = y;
	out.z = z;
	out.w = w;
#endif
}

static void Print(const TMatrix &m){
	printf(
		"%f %f %f\n"
		"%f %f %f\n"
		"%f %f %f\n"
		, m.m00, m.m01, m.m02
		, m.m10, m.m11, m.m12
		, m.m20, m.m21, m.m22
	);
}

static void Print(const CQuaternion &q){
	printf("%f %f %f %f\n", q.x, q.y, q.z, q.w);
}

static bool Float32ToHex(const TArgInfo &info){
	const TString1D &tmp = info.m_arg;
	size_t n = tmp.size();
	if(n){
		for(size_t i=0 ; i!=n ; ++i){
			float x = atof(tmp[i].c_str());
			if(0<i){
				printf(" ");
			}
			//printf("[%f]", x);
			printf("0x%08x", *(uint32*)&x);
		}
		printf("\n");
	}
	if(n==0){
		return false;
	}
	return true;
}

static bool Float64ToHex(const TArgInfo &info){
	const TString1D &tmp = info.m_arg;
	size_t n = tmp.size();
	if(n){
		assert(sizeof(double)==sizeof(uint64));
		for(size_t i=0 ; i!=n ; ++i){
			double x = atof(tmp[i].c_str());
			if(0<i){
				printf(" ");
			}
			//printf("[%f]", x);
			printf("0x%016I64x", *(uint64*)&x);
		}
		printf("\n");
	}
	if(n==0){
		return false;
	}
	return true;
}

static bool HexToFloat32(const TArgInfo &info){
	const TString1D &tmp = info.m_arg;
	size_t n = tmp.size();
	if(n){
		for(size_t i=0 ; i!=n ; ++i){
			uint32 x = strtoul(tmp[i].c_str(),NULL,16);
			if(0<i){
				printf(" ");
			}
			printf("%f", *(float*)&x);
		}
		printf("\n");
	}
	if(n==0){
		return false;
	}
	return true;
}

static bool HexToFloat64(const TArgInfo &info){
	const TString1D &tmp = info.m_arg;
	size_t n = tmp.size();
	if(n){
		for(size_t i=0 ; i!=n ; ++i){
			uint64 x = _strtoi64(tmp[i].c_str(),NULL,16);
			if(0<i){
				printf(" ");
			}
			//printf("[%016I64x]", x);
			printf("%f", *(double*)&x);
		}
		printf("\n");
	}
	if(n==0){
		return false;
	}
	return true;
}

static bool ConvertNumber(const TArgInfo &info, int dst, int src){
	const TString1D &tmp = info.m_arg;
	size_t n = tmp.size();
	if(n){
		char	buf[256];
		for(size_t i=0 ; i!=n ; ++i){
			uint64 x = _strtoi64(tmp[i].c_str(),NULL,src);
			if(0<i){
				printf(" ");
			}
			_ui64toa_s(x,buf,sizeof(buf),dst);
			if(dst==16){
				printf("0x");
			}
			printf("%s",buf);
		}
		printf("\n");
	}
	if(n==0){
		return false;
	}
	return true;
}

static bool RadianToDegree(const TArgInfo &info){
	const TString1D &tmp = info.m_arg;
	if(! (1<=tmp.size())){
		return false;
	}

	TString1D::const_iterator first = tmp.begin();
	TString1D::const_iterator last 	= tmp.end();
	for(; first!=last ; ++first){
		printf("%f\n", rad2deg(atof(first->c_str())));
	}
	return true;
}

static bool DegreeToRadian(const TArgInfo &info){
	const TString1D &tmp = info.m_arg;
	if(! (1<=tmp.size())){
		return false;
	}
	TString1D::const_iterator first = tmp.begin();
	TString1D::const_iterator last 	= tmp.end();
	for(; first!=last ; ++first){
		printf("%f\n", deg2rad(atof(first->c_str())));
	}
	return true;
}

/*
引数：幅　高さ　水平画角
（備考）「幅　高さ」はアスペクト比として使用します。
       例１		640 480
	   例２		1.333 1
*/
static bool FovVerticalToHorizontal(const TArgInfo &info){
	const TString1D &tmp = info.m_arg;

	double w = atof(tmp[0].c_str());
	double h = atof(tmp[1].c_str());
	double fov = deg2rad(atof(tmp[2].c_str()));
	double aspect = w/h;
	double ans = rad2deg(2.0 * atan (aspect * tan(fov / 2.0) ));
	printf("%f\n",ans);
	return true;
}

/*
引数：FovVerticalToHorizontalと同じ
*/
static bool FovHorizontalToVertical(const TArgInfo &info){
	const TString1D &tmp = info.m_arg;

	double w = atof(tmp[0].c_str());
	double h = atof(tmp[1].c_str());
	double fov = deg2rad(atof(tmp[2].c_str()));
	double aspect = w/h;
	double ans = rad2deg(2.0 * atan (tan(fov / 2.0) / aspect));
	printf("%f\n", ans);
	return true;
}

/*
引数：x y z q
*/
static bool QuatLength(const TArgInfo &info){
	const TString1D &tmp = info.m_arg;
	CQuaternion q(tmp[0], tmp[1], tmp[2], tmp[3]);
	printf("%f\n",q.Length());
	return true;
}

/*
引数：x y z q
*/
static bool QuatNormalize(const TArgInfo &info){
	const TString1D &tmp = info.m_arg;
	CQuaternion q(tmp[0], tmp[1], tmp[2], tmp[3]);
	float l = q.Length();
	//if(l<EPSILON){
	//}
	float d = 1.f / l;
	q.x *= d;
	q.y *= d;
	q.z *= d;
	q.w *= d;
	Print(q);
	return true;
}

/*
引数：x y z q
*/
static bool QuatInverse(const TArgInfo &info){
	const TString1D &tmp = info.m_arg;

	CQuaternion q(tmp[0], tmp[1], tmp[2], tmp[3]);
    float d = 1.f/q.Length();
    q.x *= -d;
	q.y *= -d;
	q.z *= -d;
	q.w *= d;
	Print(q);
	return true;
}

static bool QuatToMat(const TArgInfo &info){
	const TString1D &tmp = info.m_arg;

	CQuaternion q(tmp[0], tmp[1], tmp[2], tmp[3]);
	TMatrix m;
	QuatToMatMain(m,q);
	Print(m);
	return true;
}

static bool QuatToAAng(const TArgInfo &info){
	const TString1D &tmp = info.m_arg;

	CQuaternion q(tmp[0], tmp[1], tmp[2], tmp[3]);
	TAxisAngle	a;
	QuatToAxisAngleMain(a, q);
	Print(a);
	return true;
}

static bool QuatToXYZ(const TArgInfo &info){
	const TString1D &tmp = info.m_arg;

	CQuaternion q(tmp[0], tmp[1], tmp[2], tmp[3]);
	printf("あとで\n");
	return false;
}

static bool Mat33Transpose(const TArgInfo &info){
	TMatrix		m(info.m_arg);
	m.Transpose();
	Print(m);
	return true;
}

static bool Mat33Determinant(const TArgInfo &info){
	TMatrix		m(info.m_arg);
	float d = m.Determinant();
	printf("%f\n",d);
	return true;
}

static bool Mat33Invert(const TArgInfo &info){
	TMatrix		m(info.m_arg);
	m.Invert();
	Print(m);
	return true;
}

static bool Mat33Scale(const TArgInfo &info){
	TMatrix		m(info.m_arg);
	float 		x=0.f,
				y=0.f,
				z=0.f;

	m.ExtractScale(x,y,z);
	printf("%f %f %f\n", x, y, z);
	return true;
}

/*
引数:deg_x
*/
static bool Mat33RotX(const TArgInfo &info){
	const TString1D &tmp = info.m_arg;
	if(! (1 <= tmp.size())){
		return false;
	}
	TMatrix		m;
	m.RotX(deg2rad(atof(tmp[0].c_str())));
	Print(m);
	return true;
}

/*
引数:deg_y
*/
static bool Mat33RotY(const TArgInfo &info){
	const TString1D &tmp = info.m_arg;
	if(! (1 <= tmp.size())){
		return false;
	}
	TMatrix		m;
	m.RotY(deg2rad(atof(tmp[0].c_str())));
	Print(m);
	return true;
}

/*
引数:deg_z
*/
static bool Mat33RotZ(const TArgInfo &info){
	const TString1D &tmp = info.m_arg;
	if(! (1 <= tmp.size())){
		return false;
	}
	TMatrix		m;
	m.RotZ(deg2rad(atof(tmp[0].c_str())));
	Print(m);
	return true;
}

static bool Mat33RoCommon(TMatrix&mx, TMatrix&my, TMatrix&mz, const TArgInfo &info){
	const TString1D &tmp = info.m_arg;
	if(! (3 <= tmp.size())){
		return false;
	}
	mx.RotX(deg2rad(atof(tmp[0].c_str())));
	my.RotY(deg2rad(atof(tmp[1].c_str())));
	mz.RotZ(deg2rad(atof(tmp[2].c_str())));
	return true;
}

/*
引数:deg_x deg_y deg_z
*/
static bool Mat33RotYXZ(const TArgInfo &info){
	TMatrix m, mx, my, mz;
	if(! Mat33RoCommon(mx, my, mz, info)){
		return false;
	}
	m=mz*mx*my;
	Print(m);
	return true;
}

/*
引数:deg_x deg_y deg_z
*/
static bool Mat33RotZXY(const TArgInfo &info){
	TMatrix m, mx, my, mz;
	if(! Mat33RoCommon(mx, my, mz, info)){
		return false;
	}
	m=my*mx*mz;
	Print(m);
	return true;
}

/*
引数:deg_x deg_y deg_z
*/
static bool Mat33RotZYX(const TArgInfo &info){
	TMatrix m, mx, my, mz;
	if(! Mat33RoCommon(mx, my, mz, info)){
		return false;
	}
	m=mx*my*mz;
	Print(m);
	return true;
}

/*
引数:deg_x deg_y deg_z
*/
static bool Mat33RotYZX(const TArgInfo &info){
	TMatrix m, mx, my, mz;
	if(! Mat33RoCommon(mx, my, mz, info)){
		return false;
	}
	m=mx*mz*my;
	Print(m);
	return true;
}

/*
引数:deg_x deg_y deg_z
*/
static bool Mat33RotXZY(const TArgInfo &info){
	TMatrix m, mx, my, mz;
	if(! Mat33RoCommon(mx, my, mz, info)){
		return false;
	}
	m=my*mz*mx;
	Print(m);
	return true;
}

/*
引数:deg_x deg_y deg_z
*/
static bool Mat33RotXYZ(const TArgInfo &info){
	TMatrix m, mx, my, mz;
	if(! Mat33RoCommon(mx, my, mz, info)){
		return false;
	}
	m=mz*my*mx;
	Print(m);
	return true;
}

static bool Mat33EulerCommon(const TArgInfo &info, int EulOrd, const char*rot){
	TMatrix		m(info.m_arg);
	HMatrix		r;

	ConvertMatrix(r,m);
	EulerAngles angles = Eul_FromHMatrix(r, EulOrd);
	printf("%s=%f %f %f\n", rot, rad2deg(angles.x), rad2deg(angles.y), rad2deg(angles.z));
	return true;
}

static bool Mat33EulerYXZ(const TArgInfo &info){
	return Mat33EulerCommon(info, EulOrdYXZs,"yxz");
}
static bool Mat33EulerZXY(const TArgInfo &info){
	return Mat33EulerCommon(info, EulOrdZXYs,"zxy");
}
static bool Mat33EulerZYX(const TArgInfo &info){
	return Mat33EulerCommon(info, EulOrdZYXs,"zyx");
}
static bool Mat33EulerYZX(const TArgInfo &info){
	return Mat33EulerCommon(info, EulOrdYZXs,"yzx");
}
static bool Mat33EulerXZY(const TArgInfo &info){
	return Mat33EulerCommon(info, EulOrdXZYs,"xzy");
}
static bool Mat33EulerXYZ(const TArgInfo &info){
	return Mat33EulerCommon(info, EulOrdXYZs,"xyz");
}

static void Mat33ToQuat(const TArgInfo &info){
	TMatrix	m(info.m_arg);
	CQuaternion q;
	Mat33ToQuatMain(q,m);
	printf("xyxw="); Print(q);
}

static void Mat33ToAAng(const TArgInfo &info){
	const TString1D &tmp = info.m_arg;
	TMatrix		m(tmp);
	TAxisAngle	a;
	Mat33ToAAngMain(a,m);
	printf("xyxw="); Print(a);
}

static void ParseArg(TArgInfo &out, int argc, _TCHAR* argv[]){
	out.m_mode.clear();
	out.m_arg.clear();

	if(argc <= 2){
		return ;
	}

	out.m_mode.assign(argv[1]);
	for(int i=2 ; i<argc ; ++i){
		out.m_arg.push_back(argv[i]);
	}
}

int _tmain(int argc, _TCHAR* argv[])
{
	bool		result = false;
	TArgInfo	info;

	ParseArg(info,argc,argv);
	const std::string &m = info.m_mode;

	//printf("mode=%s\n", m.c_str());

	if("f32_hex" == m){
		result = Float32ToHex(info);
	}else if("f64_hex" == m){
		result = Float64ToHex(info);
	}else if("hex_f32" == m){
		result = HexToFloat32(info);
	}else if("hex_f64" == m){
		result = HexToFloat64(info);
	}else if("2_10" == m){
		result = ConvertNumber(info,10,2);
	}else if("2_16" == m){
		result = ConvertNumber(info,16,2);
	}else if("10_2" == m){
		result = ConvertNumber(info,2,10);
	}else if("10_16" == m){
		result = ConvertNumber(info,16,10);
	}else if("16_2" == m){
		result = ConvertNumber(info,2,16);
	}else if("16_10" == m){
		result = ConvertNumber(info,10,16);
	}else if("r_d" == m){
		result = RadianToDegree(info);
	}else if("d_r" == m){
		result = DegreeToRadian(info);
	}else if("v_h" == m){
		result = MatrixSub(info,FovVerticalToHorizontal,3);
	}else if("h_v" == m){
		result = MatrixSub(info,FovHorizontalToVertical,3);
	}else if("ql" == m){
		result = MatrixSub(info,QuatLength,4);
	}else if("qn" == m){
		result = MatrixSub(info,QuatNormalize,4);
	}else if("qi" == m){
		result = MatrixSub(info,QuatInverse,4);
	}else if("q_m33" == m){
		result = MatrixSub(info,QuatToMat,4);
	}else if("q_a" == m){
		result = MatrixSub(info,QuatToAAng,4);
	}else if("q_xyz" == m){
		result = MatrixSub(info,QuatToXYZ,4);
	}else if("m33_q" == m){
		result = MatrixSub(info,Mat33ToQuat,9);
	}else if("m33_a" == m){
		result = MatrixSub(info,Mat33ToAAng,9);
	}else if("m33t" == m){
		result = MatrixSub(info,Mat33Transpose,9);
	}else if("m33d" == m){
		result = MatrixSub(info,Mat33Determinant,9);
	}else if("m33i" == m){
		result = MatrixSub(info,Mat33Invert,9);
	}else if("m33_s" == m){
		result = MatrixSub(info,Mat33Scale,9);
	}else if("x_m33" == m){
		result = Mat33RotX(info);
	}else if("y_m33" == m){
		result = Mat33RotY(info);
	}else if("z_m33" == m){
		result = Mat33RotZ(info);
	}else if("yxz_m33" == m){
		result = Mat33RotYXZ(info);
	}else if("zxy_m33" == m){
		result = Mat33RotZXY(info);
	}else if("zyx_m33" == m){
		result = Mat33RotZYX(info);
	}else if("yzx_m33" == m){
		result = Mat33RotYZX(info);
	}else if("xzy_m33" == m){
		result = Mat33RotXZY(info);
	}else if("xyz_m33" == m){
		result = Mat33RotXYZ(info);
	}else if("m33_yxz" == m){
		result = MatrixSub(info,Mat33EulerYXZ,9);
	}else if("m33_zxy" == m){
		result = MatrixSub(info,Mat33EulerZXY,9);
	}else if("m33_zyx" == m){
		result = MatrixSub(info,Mat33EulerZYX,9);
	}else if("m33_yzx" == m){
		result = MatrixSub(info,Mat33EulerYZX,9);
	}else if("m33_xzy" == m){
		result = MatrixSub(info,Mat33EulerXZY,9);
	}else if("m33_xyz" == m){
		result = MatrixSub(info,Mat33EulerXYZ,9);
	}

	if(result){
		return 0;
	}
	usage();

	return 1;
}


#if 0
//==============================================================//
 //
 //    二本のベクトルで回転出来るクォータニオンを作る
 //
 //         v0    入力、 正規化すること
 //         v1    入力、 正規化すること
 //
 //==============================================================//
 void QUATERNION::quaternion_v_set( Vec v0, Vec v1 )
 {
    Vec    out;
    f32        inner, sr;

    // cross つまり外積
    out[0] = v0[1] * v1[2] - v0[2] * v1[1];
    out[1] = v0[2] * v1[0] - v0[0] * v1[2];
    out[2] = v0[0] * v1[1] - v0[1] * v1[0];
    // dot  つまり内積
    inner = v0[0] * v1[0] + v0[1] * v1[1] + v0[2] * v1[2];

    sr = sqrtf(( 1.0f + inner ) * 2.0f );
    q[0] = out[0] / sr;
    q[1] = out[1] / sr;
    q[2] = out[2] / sr;
    q[3] = sr * 0.5f;

    return;
 }
#endif

#if 0
回転前の座標が x1, y1, z1
回転後の座標が x2, y2, z2 だったとすると

double Angle = acos((x1 * x2 + y1 * y2 + z1 * z2) / ((x1 * x1 + y1 * y1 + z1 * z1) * (x2 * x2 + y2 * y2 + z2 * z2)));

でAngleの中に角度がラジアン単位で入る。

回転軸ベクトルは、回転前の位置、回転後の位置、原点（回転の中心点）の3点から、普通にポリゴンの法線ベクトル求める方法で垂直な方向ベクトルを外積で計算してやるだけ。
#endif
