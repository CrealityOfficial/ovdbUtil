#include "ovdbutil/init.h"
#include "openvdb/openvdb.h"

namespace ovdbutil
{
	void init()
	{
		openvdb::v9_0::initialize();
	}
}
