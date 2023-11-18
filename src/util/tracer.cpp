#include "tracer.h"

namespace ovdbutil
{
    TracerInterrupter::TracerInterrupter(ccglobal::Tracer* tracer)
        :m_tracer(tracer)
    {

    }

    TracerInterrupter::~TracerInterrupter()
    {

    }

    bool TracerInterrupter::wasInterrupted(int percent)
    {
        if (m_tracer)
        {
            m_tracer->progress((percent >= 0) ? ((float)percent / 100.0f) : 0.0f);
            return m_tracer->interrupt();
        }

        return false;
    }
}