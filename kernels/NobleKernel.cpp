/******************************************************************************
 * Copyright (c) 2017, Bradley J Chambers (brad.chambers@gmail.com)
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

#include "NobleKernel.hpp"

#include <memory>

#include <filters/MergeFilter.hpp>
#include <io/BufferReader.hpp>
#include <pdal/KDIndex.hpp>
#include <pdal/PDALUtils.hpp>
#include <pdal/PointView.hpp>
#include <pdal/pdal_config.hpp>
#include <pdal/pdal_macros.hpp>

namespace pdal
{

static PluginInfo const s_info = PluginInfo("kernels.noble", "Noble Kernel",
                                            "http://pdal.io/apps/noble.html");

CREATE_STATIC_PLUGIN(1, 0, NobleKernel, Kernel, s_info)

std::string NobleKernel::getName() const
{
    return s_info.name;
}

void NobleKernel::addSwitches(ProgramArgs& args)
{
    Arg& layerA = args.add("layerA", "Layer A filename", m_layerA);
    layerA.setPositional();
    Arg& layerB = args.add("layerB", "Layer B filename", m_layerB);
    layerB.setPositional();
    Arg& output = args.add("output", "Output filename", m_outputFile);
    output.setPositional();
    Arg& tolerance = args.add("tolerance", "Tolerance", m_tolerance);
    tolerance.setPositional();
}

PointViewPtr NobleKernel::loadSet(const std::string& filename,
                                  PointTable& table)
{
    Stage& reader = makeReader(filename, "");
    reader.prepare(table);
    PointViewSet viewSet = reader.execute(table);
    assert(viewSet.size() == 1);
    return *viewSet.begin();
}

int NobleKernel::execute()
{
    PointTable tableA;
    PointViewPtr viewA = loadSet(m_layerA, tableA);

    PointTable tableB;
    PointViewPtr viewB = loadSet(m_layerB, tableB);

    PointViewPtr outViewA = viewA->makeNew();
    PointViewPtr outViewB = viewB->makeNew();

    KD3Index indexA(*viewA);
    indexA.build();

    KD3Index indexB(*viewB);
    indexB.build();

    for (PointId i = 0; i < viewA->size(); ++i)
    {
        std::vector<PointId> indices(1);
        std::vector<double> sqr_dists(1);
        PointRef pointA = viewA->point(i);
        indexB.knnSearch(pointA, 1, &indices, &sqr_dists);

        if (sqr_dists[0] > m_tolerance)
            outViewA->appendPoint(*viewA, i);
    }

    for (PointId i = 0; i < viewB->size(); ++i)
    {
        std::vector<PointId> indices(1);
        std::vector<double> sqr_dists(1);
        PointRef pointB = viewB->point(i);
        indexA.knnSearch(pointB, 1, &indices, &sqr_dists);

        if (sqr_dists[0] > m_tolerance)
            outViewB->appendPoint(*viewB, i);
    }

    BufferReader bufferReader;
    bufferReader.addView(outViewA);
    bufferReader.addView(outViewB);

    // PointTable table;
    Stage& writer = makeWriter(m_outputFile, bufferReader, "");
    writer.prepare(tableA);
    writer.execute(tableA);
    return 0;
}

} // namespace pdal
