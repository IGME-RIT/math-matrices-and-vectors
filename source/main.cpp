/*
Title: Matrix Mathematics
File Name: main.cpp
Copyright © 2016
Author: Andrew Litfin
Written under the supervision of David I. Schwartz, Ph.D., and
supported by a professional development seed grant from the B. Thomas
Golisano College of Computing & Information Sciences
(https://www.rit.edu/gccis) at the Rochester Institute of Technology.
This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or (at
your option) any later version.
This program is distributed in the hope that it will be useful, but
WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
General Public License for more details.
You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

// The primary objects of study in linear algebra are matrices.
// This tutorial series will explore the applications of matrices to computer games and simulation,
//  especially in the realm of physical transformations.
// The exposition follows that of Eric Lengyel in "Foundations of Game Engine Development" (Volume 1).
// We have included the Vector structs from the previous series, and introduced Matrix structs that act similarly.
// These structs are based upon and largely follow code samples given in FGED.
//  As before, Matrix2D is heavily annotated, with other structs being annotated in places of difference.

// This tutorial explains the relationship between matrices and vectors.
// Through this, we will learn why we care about matrices at all.
// We will also revisit the definition given at the start of the last tutorial, and see why it is true.

// This tutorial is almost entirely theoretical, and can be skipped, but the author recommends against this,
//  as a full and complete knowledge of linear algebra can only make you a better programmer.

#include "../header/Matrix4D.h"
#include "../header/tests.h"
#include "../header/helpers.h"

#include <iostream>
#include <ctime>

int main()
{
	// Required for some helper functions
	srand((unsigned)time(0));

	// Let's step back from matrices for a moment and consider an interesting problem:
	// Say I have a 2D vector and I want to double its x-coordinate, but leave the y-coordinate untouched.
	// In terms of a function f : R2 -> R2 (function name f, domain R2, codomain R2), we want
	//  f(x, y) = (2x, y)
	// We can represent this as a matrix, namely
	// [ 2 0 ]
	// [ 0 1 ]
	// Now to use this matrix, we use matrix multiplication when considering the vector as a column vector.
	// Then
	//         | [ x ]
	//         | [ y ]
	//---------+-------
	// [ 2 0 ] | [ 2x ]
	// [ 0 1 ] | [ y ]
	// Hussah!

	Matrix2D m(2, 0, 0, 1);
	Vector2D x(randIntF(-10, 10), randIntF(-10, 10));
	std::cout << "x = " << x << "\n"
		"m =\n" << m << "m * x = " << m * x << "\n";

	// Now let's revisit the original definition:
	// A matrix is a way to represent any linear map between two modules,
	//  including finite-dimensional vector spaces with a basis defined for each space.
	//  In our case, these vector spaces are almost always Euclidean spaces, particularly R2, R3, or R4.
	//  Then a m by n matrix defines a function f from Rn to Rm by f(x)=Ax satisfying
	//  f(ax_1 + x_2) = af(x_1) + f(x_2), i.e. is additive and homogeneous of degree 1.

	// Ignoring the first part (take a course on linear algebra to understand), let's focus on the following:
	//  Then a m by n matrix defines a function f from Rn to Rm by f(x)=Ax satisfying
	//  f(ax_1 + x_2) = af(x_1) + f(x_2)

	// For example, a 2x2 matrix represents a function from R2 to R2.
	// But what about that other bit?
	// Say we have a 2x2 matrix A = [ [ a, b ], [ c, d ] ] and a 2D vector x=(x, y).
	// Then A*x is the vector (ax + by, cx + dy).
	// Now suppose we have another 2D vector u=(u, v).
	// Then A*u is (au + bv, cu + dv).
	// Now suppose we have a constant s.
	// Then A*(sx) = (asx + bsy, csx + dsy) = s(A*x)
	// Lastly, A*(sx + u) = (a(sx + u) + b(sy + v), c(sx + u) + d(sy + v)) = (asx + bsy, csx + dsy) + (au + bv, cu + dv) = s(A*x) + A*u.
	// This argument can be extended to an arbitrary number of dimensions.

	Matrix2D A(randIntF(-10, 10), randIntF(-10, 10), randIntF(-10, 10), randIntF(-10, 10));
	Vector2D u(randIntF(-10, 10), randIntF(-10, 10));
	float s = randFloat(-10, 10);
	std::cout << "A =\n" << A;
	std::cout << "u = " << u << ", s = " << s << "\n";
	std::cout << "A*(sx + u) = " << A*(s*x + u) << "\n";
	std::cout << "s*(A*u) + A*v = " << s*(A*x) + A*u << "\n";

	// So what does it *mean*?
	// What it means is that *the image of a vector space under a matrix is itself a vector space*.
	// That is, you can treat vectors both before and after multiplication by a matrix as elements of a vector space,
	//  and treat them as such, with the matrix providing a natural transformation between the two.

	// Corollary:
	//  The product of a matrix and a column vector is a linear combination of the columns of the matrix.
	// Proof:
	//  For any vector v in Rn with the standard basis { e_1, e_2, ... e_n } and an n by n matrix A, the product A*v can be decomposed:
	//  A*v = A*(v_1*e_1 + v_2*e_2 + ... + v_n*e_n)
	//      = v_1(A*e_1) + v_2(A*e_2) + ... + v_n(A*e_n)
	//      = v_1*A[1] + v_2*A[2] + ... + v_n*A[n] (where here we are referencing the columns by 1-indexed notation)
	
	// (Note that because of how we index matrices in code, the last line translates to code as A[0], A[1], ..., A[n-1].)

	// Indeed, from the corollary, this is how we use matrix multiplcation in code.
	// It is equivalent and smaller than the more verbose
	// Vector3D operator*(Matrix3D M, Vector3D v)
	// {
	//   return Vector3D(M(0, 0) * v.x + M(0, 1) * v.y + M(0, 2) * v.z,
	//                   M(1, 0) * v.x + M(1, 1) * v.y + M(1, 2) * v.z,
	//                   M(2, 0) * v.x + M(2, 1) * v.y + M(2, 2) * v.z);
	// }
	// Some say this means there is less chance for error, and I am inclined to agree.
	// To demonstrate, on the left side we have the explicit expansion of matrix multiplication,
	// and on the right is the linear-combination-of-columns approach.
	if (Vector2D(A(0, 0) * x.x + A(0, 1) * x.y, A(1, 0) * x.x + A(1, 1) * x.y) == (x.x * A[0] + x.y * A[1]))
	{
		std::cout << "The formulations are equivalent!\n";
	}

	// Now it is possible to multiply a row vector times a matrix, and it is important to note that the two are NOT equivalent.
	std::cout << "A*x = " << A * x << "\n";
	std::cout << "x*A = " << x * A << "\n";

	// They are equal, however, equal if you transpose the matrix, i.e.
	std::cout << "A*x (as a column vector) = x*A^T (as a row vector) = " << x * Transpose(A) << "\n";

	// This is because for all matrices A and B, if A*B is defined, then (AB)^T = B^T A^T

	std::cout << "\nPress Enter to exit . . . ";
	std::cin.get();
	return 0;
}
