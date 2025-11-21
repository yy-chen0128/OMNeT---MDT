/*
 * MDTData.cc
 *
 *  Created on: Aug 10, 2025
 *      Author: yychen
 */

#include "../../mysrc/routing/MDTData.h"

namespace mysrc {
namespace routing {

std::string MDTData::str() const
{
    std::ostringstream out;
    out << "isActive = " << isActive()
        << ", hasValidDestNum = " << hasValidDestNum()
        << ", destNum = " << getDestSeqNum()
        << ", lifetime = " << getLifeTime();

    //MDT

    const std::set<L3Address>& preList = getPrecursorList();

    if (!preList.empty()) {
        out << ", precursor list: ";
        std::set<L3Address>::const_iterator iter = preList.begin();
        out << *iter;
        for (++iter; iter != preList.end(); ++iter)
            out << "; " << *iter;
    }
    return out.str();
};

} //
} //


