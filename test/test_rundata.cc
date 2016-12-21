//
// Tests of the run data class

#include "gtest/gtest.h"

#include "dataset.h"
#include "setup.h"
#include "easylog.h"

#include <fstream>

namespace
{

// The fixture for testing class Foo.
class RunDataTest: public ::testing::Test
{
protected:
	RunDataTest()
	{
	}

	virtual ~RunDataTest()
	{
	}

	virtual void SetUp()
	{
	}

	virtual void TearDown()
	{
	}

	bool FloatEq(double d1, double d2, double epsilon = 0.001)
	{
		double diff = d1 - d2;
		if (diff < 0)
			diff = -diff;
		return diff < epsilon;
	}
};

// Tests concatenated data sets
TEST_F(RunDataTest, ConcatenatedData)
{
	int NTIMES = 10; // needs to be even
	int VSIZE = 5;
	float VAL = 7.32;

	// Create coordinates and data matrices
	NEWMAT::Matrix voxelCoords, data1, data2, data3;
	data1.ReSize(NTIMES, VSIZE * VSIZE * VSIZE);
	data2.ReSize(NTIMES, VSIZE * VSIZE * VSIZE);
	data3.ReSize(NTIMES, VSIZE * VSIZE * VSIZE);
	voxelCoords.ReSize(3, VSIZE * VSIZE * VSIZE);
	int v = 1;
	for (int z = 0; z < VSIZE; z++)
	{
		for (int y = 0; y < VSIZE; y++)
		{
			for (int x = 0; x < VSIZE; x++)
			{
				voxelCoords(1, v) = x;
				voxelCoords(2, v) = y;
				voxelCoords(3, v) = z;
				for (int n = 0; n < NTIMES; n++)
				{
					data1(n + 1, v) = VAL;
					data2(n + 1, v) = VAL * 2;
					data3(n + 1, v) = VAL * 3;
				}
				v++;
			}
		}
	}

	FabberIoMemory io;
	FabberRunData rundata(&io);
	io.SetVoxelCoords(voxelCoords);
	io.SetVoxelData("data1", data1);
	io.SetVoxelData("data2", data2);
	io.SetVoxelData("data3", data3);
	rundata.Set("data-order", "concatenate");
	NEWMAT::Matrix data = rundata.GetMainVoxelData();

	ASSERT_EQ(data.Nrows(), NTIMES * 3);
	ASSERT_EQ(data.Ncols(), VSIZE * VSIZE * VSIZE);
	for (int i = 0; i < VSIZE * VSIZE * VSIZE; i++)
	{
		for (int t = 0; t < NTIMES * 3; t++)
		{
			if (t < NTIMES)
				ASSERT_FLOAT_EQ(data(t + 1, i + 1), VAL);
			else if (t < NTIMES * 2)
				ASSERT_FLOAT_EQ(data(t + 1, i + 1), VAL * 2);
			else
				ASSERT_FLOAT_EQ(data(t + 1, i + 1), VAL * 3);
		}
	}
}

// Tests interleaved data sets
TEST_F(RunDataTest, InterleavedData)
{
	int NTIMES = 10; // needs to be even
	int VSIZE = 5;
	float VAL = 7.32;

	// Create coordinates and data matrices
	NEWMAT::Matrix voxelCoords, data1, data2, data3;
	data1.ReSize(NTIMES, VSIZE * VSIZE * VSIZE);
	data2.ReSize(NTIMES, VSIZE * VSIZE * VSIZE);
	data3.ReSize(NTIMES, VSIZE * VSIZE * VSIZE);
	voxelCoords.ReSize(3, VSIZE * VSIZE * VSIZE);
	int v = 1;
	for (int z = 0; z < VSIZE; z++)
	{
		for (int y = 0; y < VSIZE; y++)
		{
			for (int x = 0; x < VSIZE; x++)
			{
				voxelCoords(1, v) = x;
				voxelCoords(2, v) = y;
				voxelCoords(3, v) = z;
				for (int n = 0; n < NTIMES; n++)
				{
					data1(n + 1, v) = VAL;
					data2(n + 1, v) = VAL * 2;
					data3(n + 1, v) = VAL * 3;
				}
				v++;
			}
		}
	}

	FabberIoMemory io;
	FabberRunData rundata(&io);
	io.SetVoxelCoords(voxelCoords);
	io.SetVoxelData("data1", data1);
	io.SetVoxelData("data2", data2);
	io.SetVoxelData("data3", data3);
	rundata.Set("data-order", "interleave");
	NEWMAT::Matrix data = rundata.GetMainVoxelData();

	ASSERT_EQ(data.Nrows(), NTIMES * 3);
	ASSERT_EQ(data.Ncols(), VSIZE * VSIZE * VSIZE);
	for (int i = 0; i < VSIZE * VSIZE * VSIZE; i++)
	{
		for (int t = 0; t < NTIMES * 3; t++)
		{
			if (t % 3 == 0)
				ASSERT_FLOAT_EQ(data(t + 1, i + 1), VAL);
			else if (t % 3 == 1)
				ASSERT_FLOAT_EQ(data(t + 1, i + 1), VAL * 2);
			else
				ASSERT_FLOAT_EQ(data(t + 1, i + 1), VAL * 3);
		}
	}
}
// Tests inconsistent multi-data sets
TEST_F(RunDataTest, MultiDataInconsistent)
{
	int NTIMES = 10; // needs to be even
	int VSIZE = 5;
	float VAL = 7.32;

	// Create coordinates and data matrices
	NEWMAT::Matrix voxelCoords, data1, data2, data3;
	data1.ReSize(NTIMES, VSIZE * VSIZE * VSIZE);
	data2.ReSize(NTIMES, VSIZE * VSIZE * VSIZE);
	data3.ReSize(NTIMES, VSIZE * VSIZE * VSIZE);
	voxelCoords.ReSize(3, VSIZE * VSIZE * VSIZE);
	int v = 1;
	for (int z = 0; z < VSIZE; z++)
	{
		for (int y = 0; y < VSIZE; y++)
		{
			for (int x = 0; x < VSIZE; x++)
			{
				voxelCoords(1, v) = x;
				voxelCoords(2, v) = y;
				voxelCoords(3, v) = z;
				for (int n = 0; n < NTIMES; n++)
				{
					data1(n + 1, v) = VAL;
					data2(n + 1, v) = VAL * 2;
					data3(n + 1, v) = VAL * 3;
				}
				v++;
			}
		}
	}

	FabberIoMemory io;
	FabberRunData rundata(&io);
	io.SetVoxelCoords(voxelCoords);
	io.SetVoxelData("data1", data1);
	io.SetVoxelData("data2", data2);
	io.SetVoxelData("data3", data3);
	rundata.Set("data-order", "singlefile");
	ASSERT_THROW(NEWMAT::Matrix data = rundata.GetMainVoxelData(), Invalid_option);
}

// Tests that we can read config from a .fab file
TEST_F(RunDataTest, OptionsFile)
{
	string FILENAME = "test_config";

	ofstream os;
	os.open(FILENAME.c_str(), ios::out);
	os << "noise=white" << endl << "model=trivial" << endl << "method=vb" << endl << "bool-option" << endl
			<< "#comment, ignored" << endl;
	os.close();

	FabberIoMemory io;
	FabberRunData rundata(&io);
	rundata.ParseParamFile(FILENAME);
	ASSERT_EQ("white", rundata.GetString("noise"));
	ASSERT_EQ("trivial", rundata.GetString("model"));
	ASSERT_EQ("vb", rundata.GetString("method"));
	ASSERT_EQ(true, rundata.GetBool("bool-option"));
}

// Tests embedded comments
TEST_F(RunDataTest, OptionsFileEmbeddedComment)
{
	string FILENAME = "test_config";

	ofstream os;
	os.open(FILENAME.c_str(), ios::out);
	os << "model=trivial #just keep things simple" << endl;
	os.close();

	FabberIoMemory io;
	FabberRunData rundata(&io);
	rundata.ParseParamFile(FILENAME);
	ASSERT_EQ("trivial", rundata.GetString("model"));
}

// Tests unsetting an option
TEST_F(RunDataTest, Unset)
{
	FabberIoMemory io;
	FabberRunData rundata(&io);
	rundata.Set("wibble", "wobble");
	rundata.SetBool("bobble");

	ASSERT_EQ("wobble", rundata.GetStringDefault("wibble", "squabble"));
	rundata.Unset("wibble");
	ASSERT_EQ("squabble", rundata.GetStringDefault("wibble", "squabble"));
	ASSERT_EQ(true, rundata.GetBool("bobble"));
	rundata.Unset("bobble");
	ASSERT_EQ(false, rundata.GetBool("bobble"));
}

// Tests odd case where data name could lead to circular reference
TEST_F(RunDataTest, CircularDataRef)
{
	int NTIMES = 10; // needs to be even
	int VSIZE = 5;
	float VAL = 7.32;

	// Create coordinates and data matrices
	NEWMAT::Matrix voxelCoords, data1;
	data1.ReSize(NTIMES, VSIZE * VSIZE * VSIZE);
	voxelCoords.ReSize(3, VSIZE * VSIZE * VSIZE);
	int v = 1;
	for (int z = 0; z < VSIZE; z++)
	{
		for (int y = 0; y < VSIZE; y++)
		{
			for (int x = 0; x < VSIZE; x++)
			{
				voxelCoords(1, v) = x;
				voxelCoords(2, v) = y;
				voxelCoords(3, v) = z;
				for (int n = 0; n < NTIMES; n++)
				{
					data1(n + 1, v) = VAL;
				}
				v++;
			}
		}
	}
	FabberIoMemory io;
	io.SetVoxelCoords(voxelCoords);
	io.SetVoxelData("data", data1);

	FabberRunData rundata(&io);
	rundata.Set("data", "data");

	NEWMAT::Matrix data = rundata.GetMainVoxelData();

	ASSERT_EQ(VSIZE * VSIZE * VSIZE, data.Ncols());
	ASSERT_EQ(NTIMES, data.Nrows());
}

// Tests clearing voxel data
TEST_F(RunDataTest, ClearVoxelData)
{
	int NTIMES = 10; // needs to be even
	int VSIZE = 5;
	float VAL = 7.32;

	// Create coordinates and data matrices
	NEWMAT::Matrix voxelCoords, data1, data2, data3;
	data1.ReSize(NTIMES, VSIZE * VSIZE * VSIZE);
	data2.ReSize(NTIMES, VSIZE * VSIZE * VSIZE);
	data3.ReSize(NTIMES, VSIZE * VSIZE * VSIZE);
	voxelCoords.ReSize(3, VSIZE * VSIZE * VSIZE);
	int v = 1;
	for (int z = 0; z < VSIZE; z++)
	{
		for (int y = 0; y < VSIZE; y++)
		{
			for (int x = 0; x < VSIZE; x++)
			{
				voxelCoords(1, v) = x;
				voxelCoords(2, v) = y;
				voxelCoords(3, v) = z;
				for (int n = 0; n < NTIMES; n++)
				{
					data1(n + 1, v) = VAL;
					data2(n + 1, v) = VAL * 2;
					data3(n + 1, v) = VAL * 3;
				}
				v++;
			}
		}
	}

	FabberIoMemory io;
	FabberRunData rundata(&io);
	io.SetVoxelCoords(voxelCoords);
	io.SetVoxelData("data1", data1);
	io.SetVoxelData("data2", data2);
	io.SetVoxelData("data3", data3);

	io.ClearVoxelData("data1");
	ASSERT_NO_THROW(rundata.GetVoxelCoords());
	ASSERT_THROW(rundata.GetVoxelData("data1"), DataNotFound);
	ASSERT_NO_THROW(rundata.GetVoxelData("data2"));
	ASSERT_NO_THROW(rundata.GetVoxelData("data3"));

	io.ClearVoxelData();
	ASSERT_NO_THROW(rundata.GetVoxelCoords());
	ASSERT_THROW(rundata.GetVoxelData("data1"), DataNotFound);
	ASSERT_THROW(rundata.GetVoxelData("data2"), DataNotFound);
	ASSERT_THROW(rundata.GetVoxelData("data3"), DataNotFound);
}

}
