/** This is a test application for GeoWave Math library. */

#include "../GW_Maths.h"
#include "../GW_Vector2D.h"
#include "../GW_Vector3D.h"
#include "../GW_Vector4D.h"
#include "../GW_VectorND.h"
#include "../GW_Matrix2x2.h"
#include "../GW_Matrix3x3.h"
#include "../GW_Matrix4x4.h"
#include "../GW_MatrixNxP.h"
#include "../GW_SparseMatrix.h"
#include "../GW_MatrixStatic.h"
#include "../GW_Quaternion.h"

using namespace GW;

void main()
{
	cout << "############################################################" << endl;
	cout << "          GeoWave Maths Library  --  Tests" << endl;
	cout << "############################################################" << endl;
	GW_Maths::TestClass();
	GW_VectorStatic<3,GW_Float>::TestClass();
	GW_Vector2D::TestClass();
	GW_Vector3D::TestClass();
	GW_Vector4D::TestClass();
	GW_VectorND::TestClass();
	GW_MatrixStatic<3,3,GW_Float>::TestClass();
	GW_Matrix2x2::TestClass();
	GW_Matrix3x3::TestClass();
	GW_Matrix4x4::TestClass();
	GW_MatrixNxP::TestClass();
	GW_SparseMatrix::TestClass();
	GW_Quaternion::TestClass();
	cout << "############################################################" << endl;
}