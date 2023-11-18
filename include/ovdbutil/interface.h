#ifndef OVDBUTIL_INTERFACE_1604911737496_H
#define OVDBUTIL_INTERFACE_1604911737496_H
#include "ccglobal/export.h"

#if USE_OVDBUTIL_DLL
	#define OVDBUTIL_API CC_DECLARE_IMPORT
#elif USE_OVDBUTIL_STATIC
	#define OVDBUTIL_API CC_DECLARE_STATIC
#else
	#if OVDBUTIL_DLL
		#define OVDBUTIL_API CC_DECLARE_EXPORT
	#else
		#define OVDBUTIL_API CC_DECLARE_STATIC
	#endif
#endif

#include "trimesh2/TriMesh.h"
#include "trimesh2/XForm.h"
#include "ccglobal/tracer.h"
#include <memory>

namespace ovdbutil {
	typedef std::shared_ptr<trimesh::TriMesh> DBMeshPtr;
}
#endif // OVDBUTIL_INTERFACE_1604911737496_H