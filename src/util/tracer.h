#ifndef TRACER_OVDBUTIL_UTIL_H
#define TRACER_OVDBUTIL_UTIL_H
#include "openvdb/util/NullInterrupter.h"
#include "ccglobal/tracer.h"

namespace ovdbutil
{
    class TracerInterrupter : public openvdb::util::NullInterrupter
    {
    public:
        TracerInterrupter(ccglobal::Tracer* tracer);
        virtual ~TracerInterrupter();

        bool wasInterrupted(int /*percent*/ = -1) override;
    protected:
        ccglobal::Tracer* m_tracer;
    };
}

#endif // OVDBUTIL_UTIL_H