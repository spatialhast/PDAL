/******************************************************************************
* Copyright (c) 2016, Bradley J Chambers (brad.chambers@gmail.com)
*
* All rights reserved.
*
* Redistribution and use in source and binary forms, with or without
* modification, are permitted provided that the following
* conditions are met:
*
*     * Redistributions of source code must retain the above copyright
*       notice, this list of conditions and the following disclaimer.
*     * Redistributions in binary form must reproduce the above copyright
*       notice, this list of conditions and the following disclaimer in
*       the documentation and/or other materials provided
*       with the distribution.
*     * Neither the name of Hobu, Inc. or Flaxen Geo Consulting nor the
*       names of its contributors may be used to endorse or promote
*       products derived from this software without specific prior
*       written permission.
*
* THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
* "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
* LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS
* FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE
* COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
* INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,
* BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS
* OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED
* AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
* OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT
* OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY
* OF SUCH DAMAGE.
****************************************************************************/

#include "RadialDensityFilter.hpp"

#include <pdal/KDIndex.hpp>
#include <pdal/pdal_macros.hpp>

#include <string>
#include <vector>

namespace pdal
{

static PluginInfo const s_info =
    PluginInfo("filters.radialdensity", "RadialDensity Filter",
               "http://pdal.io/stages/filters.radialdensity.html");

CREATE_STATIC_PLUGIN(1, 0, RadialDensityFilter, Filter, s_info)

std::string RadialDensityFilter::getName() const
{
    return s_info.name;
}

void RadialDensityFilter::addArgs(ProgramArgs& args)
{
    args.add("radius", "Radius", m_rad, 1.0);
}

void RadialDensityFilter::addDimensions(PointLayoutPtr layout)
{
    using namespace Dimension;
    m_rdens = layout->registerOrAssignDim("RadialDensity", Type::Double);
}

void RadialDensityFilter::filter(PointView& view)
{
    using namespace Dimension;
    
    // Build the 3D KD-tree.
    log()->get(LogLevel::Debug) << "Building 3D KD-tree...\n";
    KD3Index index(view);
    index.build();
 
    // Search for neighboring points within the specified radius. The number of
    // neighbors (which includes the query point) is normalized by the volume
    // of the search sphere and recorded as the density.
    log()->get(LogLevel::Debug) << "Computing densities...\n";
    double factor = 1.0 / ((4.0 / 3.0) * 3.14159 * (m_rad * m_rad * m_rad));
    for (PointId i = 0; i < view.size(); ++i)
    {
        std::vector<PointId> pts = index.radius(i, m_rad);
        view.setField(m_rdens, i, pts.size() * factor);
    } 
}

} // namespace pdal
