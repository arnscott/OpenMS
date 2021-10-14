// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2018.
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
// $Maintainer: Hannes Roest $
// $Authors: Hannes Roest $
// --------------------------------------------------------------------------

#include <OpenMS/ANALYSIS/OPENSWATH/DATAACCESS/SpectrumAccessQuadMZTransforming.h>

#include <iostream>

namespace OpenMS
{

  SpectrumAccessQuadMZTransforming::SpectrumAccessQuadMZTransforming(
      OpenSwath::SpectrumAccessPtr sptr,
      double a, double b, double c, bool ppm) :
        SpectrumAccessTransforming(sptr), 
        a_(a), 
        b_(b), 
        c_(c), 
        ppm_(ppm)
    {}
        
    SpectrumAccessQuadMZTransforming::~SpectrumAccessQuadMZTransforming() {}

    boost::shared_ptr<OpenSwath::ISpectrumAccess> SpectrumAccessQuadMZTransforming::lightClone() const
    {
      // Create a light clone of *this by initializing a new
      // SpectrumAccessQuadMZTransforming with a light clone of the underlying
      // SpectrumAccess object and the parameters.
      return boost::shared_ptr<SpectrumAccessQuadMZTransforming>(
          new SpectrumAccessQuadMZTransforming(sptr_->lightClone(), a_, b_, c_, ppm_));
    }

    OpenSwath::SpectrumPtr SpectrumAccessQuadMZTransforming::getSpectrumById(int id)
    {
      static int ctr = 0;
      OpenSwath::SpectrumPtr s = sptr_->getSpectrumById(id);
      // with: OpenSwathWorkflow took 05:58 m (wall), 28:40 m (CPU), 36.62 s (system), 28:04 m (user).
      // w/o: OpenSwathWorkflow took 03:16 m (wall), 13:50 m (CPU), 21.27 s (system), 13:29 m (user).
      // with new it :OpenSwathWorkflow took 03:02 m (wall), 14:41 m (CPU), 20.05 s (system), 14:21 m (user).
      // with no calibration at all: OpenSwathWorkflow took 02:51 m (wall), 14:01 m (CPU), 18.36 s (system), 13:43 m (user).
      // now takes: OpenSwathWorkflow took 05:05 m (wall), 22:20 m (CPU), 31.41 s (system), 21:48 m (user).
      //
      // using regression_delta_ppm 3 threads: 
      //    -- done [took 08:07 m (CPU), 02:49 m (Wall)] --
      //    OpenSwathWorkflow took 03:52 m (wall), 09:47 m (CPU), 11.67 s (system), 09:35 m (user).
      // using none 3 threads: 
      //    -- done [took 08:39 m (CPU), 03:01 m (Wall)] -- 
      //    OpenSwathWorkflow took 04:25 m (wall), 11:37 m (CPU), 15.26 s (system), 11:22 m (user).
      // using old code & regression_delta_ppm 3 threads: 
      //    -- done [took 09:41 m (CPU), 03:26 m (Wall)] --
      //    OpenSwathWorkflow took 04:19 m (wall), 11:18 m (CPU), 11.05 s (system), 11:07 m (user).
      //
      // using openms v 2.4 code & 3 threads: 
      //    -- done [took 08:37 m (CPU), 03:07 m (Wall)] --
      //    OpenSwathWorkflow took 05:09 m (wall), 10:22 m (CPU), 11.27 s (system), 10:11 m (user).
      //
      // using default params: no MS1, no calibration, fewer peaks
      //
      //
      //
#if 0
      // auto it = s->getMZArray()->data.begin();
      // if (c_ == 0)
      // {
      //   for (auto& it : s->getMZArray()->data)
      //   {
      //     // mz = a + b * mz + c * mz^2
      //     double predict = 
      //       a_ + 
      //       b_ * it;

      //     // If ppm is true, we predicted the ppm deviation, not the actual new mass
      //     if (ppm_)
      //     {
      //       it = it - predict * it / 1000000;
      //     }
      //     else
      //     {
      //       it = predict;
      //     }
      //   }
      // }

      // else
      {
      for (auto& it : s->getMZArray()->data)
      {
        // mz = a + b * mz + c * mz^2
        double predict = 
          a_ + 
          b_ * it +
          c_ * it * it;

        // If ppm is true, we predicted the ppm deviation, not the actual new mass
        if (ppm_)
        {
          it = it - predict * it / 1000000;
        }
        else
        {
          it = predict;
        }
      }
      }
#else
      for (size_t i = 0; i < s->getMZArray()->data.size(); i++)
      {
        // mz = a + b * mz + c * mz^2
        double predict = 
          a_ + 
          b_ * s->getMZArray()->data[i] +
          c_ * s->getMZArray()->data[i] * s->getMZArray()->data[i];

        // If ppm is true, we predicted the ppm deviation, not the actual new mass
        if (ppm_)
        {
          s->getMZArray()->data[i] = s->getMZArray()->data[i] - predict*s->getMZArray()->data[i]/1000000;
        }
        else
        {
          s->getMZArray()->data[i] = predict;
        }
      }
#endif
      // std::cout << " corrected " << ctr++ << " spectra " << std::endl;
      //  corrected 97706 spectra  [for     <scan num="50882"
      //   -> each spectrum is corrected twice ??
      return s;
    }

}
