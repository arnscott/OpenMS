// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry               
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2021.
// 
// This software is released under a three-clause BSD license:
//  * Redistributions of source code must retain the above copyright
//    notice, this list of conditions and the following disclaimer.
//  * Redistributions in binary form must reproduce the above copyright
//    notice, this list of conditions and the following disclaimer in the
//    documentation and/or other materials provided with the distribution.
//  * Neither the name of any author or any participating institution 
//    may be used to endorse or promote products derived from this software 
//    without specific prior written permission.
// For a full list of authors, refer to the file AUTHORS. 
// --------------------------------------------------------------------------
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL ANY OF THE AUTHORS OR THE CONTRIBUTING 
// INSTITUTIONS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, 
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, 
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; 
// OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, 
// WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR 
// OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF 
// ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
// 
// --------------------------------------------------------------------------
// $Maintainer: Mathias Walzer$
// $Authors: $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

///////////////////////////
#include <OpenMS/COMPARISON/CLUSTERING/SingleLinkage.h>
#include <OpenMS/COMPARISON/CLUSTERING/ClusterAnalyzer.h>
#include <OpenMS/DATASTRUCTURES/DistanceMatrix.h>
#include <vector>
///////////////////////////

using namespace OpenMS;
using namespace std;

START_TEST(SingleLinkage, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

SingleLinkage* ptr = nullptr;
SingleLinkage* nullPointer = nullptr;
START_SECTION(SingleLinkage())
{
	ptr = new SingleLinkage();
	TEST_NOT_EQUAL(ptr, nullPointer)
}
END_SECTION

START_SECTION(~SingleLinkage())
{
	delete ptr;
}
END_SECTION

ptr = new SingleLinkage();

START_SECTION((SingleLinkage(const SingleLinkage &source)))
{
	SingleLinkage copy(*ptr);
	TEST_EQUAL(copy.getProductName(), ptr->getProductName());
}
END_SECTION

START_SECTION((SingleLinkage& operator=(const SingleLinkage &source)))
{
	SingleLinkage copy;
	copy = *ptr;
	TEST_EQUAL(copy.getProductName(), ptr->getProductName());
}
END_SECTION

START_SECTION((void operator()(DistanceMatrix< float > &original_distance, std::vector<BinaryTreeNode>& cluster_tree, const float threshold=1) const))
{

	DistanceMatrix<float> matrix(6,666);
	matrix.setValue(1,0,0.5f);
	matrix.setValue(2,0,0.8f);
	matrix.setValue(2,1,0.3f);
	matrix.setValue(3,0,0.6f);
	matrix.setValue(3,1,0.8f);
	matrix.setValue(3,2,0.8f);
	matrix.setValue(4,0,0.8f);
	matrix.setValue(4,1,0.8f);
	matrix.setValue(4,2,0.8f);
	matrix.setValue(4,3,0.4f);
	matrix.setValue(5,0,0.7f);
	matrix.setValue(5,1,0.8f);
	matrix.setValue(5,2,0.8f);
	matrix.setValue(5,3,0.8f);
	matrix.setValue(5,4,0.8f);

	vector< BinaryTreeNode > result;
	vector< BinaryTreeNode > tree;
	//~ tree.push_back(BinaryTreeNode(1,2,0.3f));
	//~ tree.push_back(BinaryTreeNode(2,3,0.4f));
	//~ tree.push_back(BinaryTreeNode(0,1,0.5f));
	//~ tree.push_back(BinaryTreeNode(0,1,0.6f));
	//~ tree.push_back(BinaryTreeNode(0,1,0.7f));
	tree.push_back(BinaryTreeNode(1,2,0.3f));
	tree.push_back(BinaryTreeNode(3,4,0.4f));
	tree.push_back(BinaryTreeNode(0,1,0.5f));
	tree.push_back(BinaryTreeNode(0,3,0.6f));
	tree.push_back(BinaryTreeNode(0,5,0.7f));

	(*ptr)(matrix,result);
	TEST_EQUAL(tree.size(), result.size());
	for (Size i = 0; i < tree.size(); ++i)
	{
			TOLERANCE_ABSOLUTE(0.0001);
			TEST_EQUAL(tree[i].left_child, result[i].left_child);
			TEST_EQUAL(tree[i].right_child, result[i].right_child);
			TEST_REAL_SIMILAR(tree[i].distance, result[i].distance);
	}

	//~ float th(0.6f);
	//~ tree.pop_back();
	//~ tree.pop_back();
	//~ result.clear();

	//~ (*ptr)(matrix,result,th);
	//~ for (Size i = 0; i < tree.size(); ++i)
	//~ {
			//~ TOLERANCE_ABSOLUTE(0.0001);
			//~ TEST_EQUAL(tree[i].left_child, result[i].left_child);
			//~ TEST_EQUAL(tree[i].right_child, result[i].right_child);
			//~ TEST_REAL_SIMILAR(tree[i].distance, result[i].distance);
	//~ }
}
END_SECTION

START_SECTION((static const String getProductName()))
{
  TEST_EQUAL(ptr->getProductName(), "SingleLinkage")
}
END_SECTION

START_SECTION((static ClusterFunctor* create()))
{
	ClusterFunctor* cf = SingleLinkage::create();
  TEST_NOT_EQUAL( dynamic_cast<SingleLinkage*>(cf) , nullPointer)
  delete cf;
}
END_SECTION

delete ptr;

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST



