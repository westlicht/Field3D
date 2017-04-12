//----------------------------------------------------------------------------//

/*
 * Copyright (c) 2009 Sony Pictures Imageworks Inc
 *
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions
 * are met:
 *
 * Redistributions of source code must retain the above copyright
 * notice, this list of conditions and the following disclaimer.
 * Redistributions in binary form must reproduce the above copyright
 * notice, this list of conditions and the following disclaimer in the
 * documentation and/or other materials provided with the
 * distribution.  Neither the name of Sony Pictures Imageworks nor the
 * names of its contributors may be used to endorse or promote
 * products derived from this software without specific prior written
 * permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 * "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
 * LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS
 * FOR A PARTICULAR PURPOSE ARE DISCLAIMED.  IN NO EVENT SHALL THE
 * COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
 * INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
 * (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
 * SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
 * HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT,
 * STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
 * ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED
 * OF THE POSSIBILITY OF SUCH DAMAGE.
 */

//----------------------------------------------------------------------------//

/*! \file Resample.h
  \brief Contains functions for resampling fields
*/

//----------------------------------------------------------------------------//

#ifndef _INCLUDED_Field3D_Resample_H_
#define _INCLUDED_Field3D_Resample_H_

#include <vector>

#include <boost/thread/thread.hpp>
#include <boost/thread/mutex.hpp>

#include "DenseField.h"
#include "SparseField.h"

//----------------------------------------------------------------------------//

/* TODO LIST

 * x Implement dumb, dense resampling
 * x For SparseField, only write non-zero results
 * x Implement more filters
 * For SparseField, be smart about which blocks are computed
 * x Multi-threading using boost
 * Multi-threading using TBB

 */

//----------------------------------------------------------------------------//

#include "ns.h"

FIELD3D_NAMESPACE_OPEN

//----------------------------------------------------------------------------//
// Resizing functions 
//----------------------------------------------------------------------------//

//! Resamples the source field into the target field, such that the 
//! new data window is @dataWindow.
//! \note This will query filter.isSeparable() and call separableResample()
//! if possible.
//! \note The extents of the field will be reset to match the data window.
//! This should 
template <typename Field_T, typename FilterOp_T>
bool resample(const Field_T &src, Field_T &tgt, const V3i &newRes,
              const FilterOp_T &filter, const size_t numThreads = 1);

//----------------------------------------------------------------------------//
// Filter
//----------------------------------------------------------------------------//

struct Filter
{
  // Typedefs ---

  typedef boost::shared_ptr<Filter>       Ptr;
  typedef boost::shared_ptr<const Filter> CPtr;

  // To be overridden by subclasses ---

  //! Evaluates the filter at coordinate 't'
  virtual float eval(const float t) const = 0;
  //! Radial width of the filter (half of diameter)
  virtual float support()           const = 0;

  // May be overridden by subclasses ---

  //! Initial value (zero by default, but need to be different for min/max)
  virtual float initialValue()      const
  { return 0.0f; }

};

//----------------------------------------------------------------------------//
// BoxFilter
//----------------------------------------------------------------------------//

struct BoxFilter : public Filter
{
  // Typedefs
  typedef boost::shared_ptr<BoxFilter>       Ptr;
  typedef boost::shared_ptr<const BoxFilter> CPtr;

  static const bool isAnalytic = false;

  // Ctors
  BoxFilter()
    : m_width(1.0)
  { }
  BoxFilter(const float width)
    : m_width(width)
  { }
  // From Filter base class 
  virtual float eval(const float x) const
  {
    const float t = x / m_width;
    if (t <= 0.5f) {
      return 1.0f;
    } else {
      return 0.0f;
    }
  }
  virtual float support() const
  { 
    return 0.5f * m_width; 
  }
  template <typename Value_T>
  static void op(Value_T &accumValue, const Value_T value) 
  { /* no-op */ }
private:
  const float m_width;
};

//----------------------------------------------------------------------------//
// MinFilter
//----------------------------------------------------------------------------//

struct MinFilter : public Filter
{
  // Typedefs
  typedef boost::shared_ptr<MinFilter>       Ptr;
  typedef boost::shared_ptr<const MinFilter> CPtr;

  static const bool isAnalytic = true;

  // Ctors
  MinFilter()
    : m_width(1.0)
  { }
  MinFilter(const float width)
    : m_width(width)
  { }
  // From Filter base class 
  virtual float eval(const float x) const
  {
    const float t = x / m_width;
    if (t <= 0.5f) {
      return 1.0f;
    } else {
      return 0.0f;
    }
  }
  virtual float support() const
  { 
    return 0.5f * m_width; 
  }
  virtual float initialValue() const
  {
    return std::numeric_limits<float>::max();
  }

  template <typename T>
  static void op(Imath::Vec3<T> &accumValue, const Imath::Vec3<T> value)
  {
    accumValue.x = std::min(accumValue.x, value.x);
    accumValue.y = std::min(accumValue.y, value.y);
    accumValue.z = std::min(accumValue.z, value.z);
  }

  template <typename Value_T>
  static void op(Value_T &accumValue, const Value_T value)
  {
    accumValue = std::min(accumValue, value);
  }

private:
  const float m_width;
};

//----------------------------------------------------------------------------//
// MaxFilter
//----------------------------------------------------------------------------//

struct MaxFilter : public Filter
{
  // Typedefs
  typedef boost::shared_ptr<MaxFilter>       Ptr;
  typedef boost::shared_ptr<const MaxFilter> CPtr;

  static const bool isAnalytic = true;

  // Ctors
  MaxFilter()
    : m_width(1.0)
  { }
  MaxFilter(const float width)
    : m_width(width)
  { }
  // From Filter base class 
  virtual float eval(const float x) const
  {
    const float t = x / m_width;
    if (t <= 0.5f) {
      return 1.0f;
    } else {
      return 0.0f;
    }
  }
  virtual float support() const
  { 
    return 0.5f * m_width; 
  }
  virtual float initialValue() const
  {
    return -std::numeric_limits<float>::max();
  }


  template <typename T>
  static void op(Imath::Vec3<T> &accumValue, const Imath::Vec3<T> value)
  {
    accumValue.x = std::max(accumValue.x, value.x);
    accumValue.y = std::max(accumValue.y, value.y);
    accumValue.z = std::max(accumValue.z, value.z);
  }

  template <typename Value_T>
  static void op(Value_T &accumValue, const Value_T value)
  {
    accumValue = std::max(accumValue, value);
  }

private:
  const float m_width;
};

//----------------------------------------------------------------------------//
// TriangleFilter
//----------------------------------------------------------------------------//

struct TriangleFilter : public Filter
{
  // Typedefs
  typedef boost::shared_ptr<TriangleFilter>       Ptr;
  typedef boost::shared_ptr<const TriangleFilter> CPtr;

  static const bool isAnalytic = false;

  // Ctors
  TriangleFilter()
    : m_width(1.0)
  { }
  TriangleFilter(const float width)
    : m_width(width)
  { }
  // From Filter base class 
  virtual float eval(const float x) const
  {
    const float t = x / m_width;
    if (t > 1.0) {
      return 0.0f;
    }
    return 1.0f - t;
  }
  virtual float support() const
  {
    return 1.0f * m_width;
  }
  template <typename Value_T>
  static void op(Value_T &/*accumValue*/, const Value_T /*value*/)
  { /* No-op */ }
private:
  const float m_width;
};

//----------------------------------------------------------------------------//
// GaussianFilter
//----------------------------------------------------------------------------//

struct GaussianFilter : public Filter
{
  // Typedefs
  typedef boost::shared_ptr<GaussianFilter>       Ptr;
  typedef boost::shared_ptr<const GaussianFilter> CPtr;

  static const bool isAnalytic = false;

  // Ctor 
  GaussianFilter(const float alpha = 2.0, const float width = 2.0)
    : m_alpha(alpha), 
      m_exp(std::exp(-alpha * width * width)),
      m_width(width)
  { /* Empty */ }
  // From Filter base class 
  virtual float eval(const float t) const
  {
    const float x = t / m_width;
    return std::max(0.0f, std::exp(-m_alpha * x * x) - m_exp);
  }
  virtual float support() const
  {
    return 2.0f * m_width;
  }
  template <typename Value_T>
  static void op(Value_T &accumValue, const Value_T value)
  { /* No-op */ }
private:
  const float m_alpha, m_exp, m_width;
};

//----------------------------------------------------------------------------//
// MitchellFilter
//----------------------------------------------------------------------------//

struct MitchellFilter : public Filter
{
  // Typedefs
  typedef boost::shared_ptr<MitchellFilter>       Ptr;
  typedef boost::shared_ptr<const MitchellFilter> CPtr;

  static const bool isAnalytic = false;

  // Ctor 
  MitchellFilter(const float width = 1.0, 
                 const float B = 1.0 / 3.0, const float C = 1.0 / 3.0)
    : m_B(B), m_C(C), m_width(width)
  { /* Empty */ }
  // From Filter base class 
  virtual float eval(const float x) const
  {
    const float ax = std::abs(x / m_width);
    if (ax < 1) {
      return ((12 - 9 * m_B - 6 * m_C) * ax * ax * ax +
              (-18 + 12 * m_B + 6 * m_C) * ax * ax + (6 - 2 * m_B)) / 6;
    } else if ((ax >= 1) && (ax < 2)) {
      return ((-m_B - 6 * m_C) * ax * ax * ax +
              (6 * m_B + 30 * m_C) * ax * ax + (-12 * m_B - 48 * m_C) *
              ax + (8 * m_B + 24 * m_C)) / 6;
    } else {
      return 0;
    }
  }
  virtual float support() const
  {
    return 2.0f * m_width;
  }
  template <typename Value_T>
  static void op(Value_T &accumValue, const Value_T value)
  { /* No-op */ }
private:
  const float m_B, m_C;
  const float m_width;
};

//----------------------------------------------------------------------------//
// Implementation details
//----------------------------------------------------------------------------//

namespace detail {

  //--------------------------------------------------------------------------//

  Box3i srcSupportBBox(const V3f &tgtP, const float support, const V3i &doUpres,
                       const V3f &srcSize, const V3f &tgtSize);

  //--------------------------------------------------------------------------//

  std::pair<int, int>
  srcSupportBBox(const float &tgtP, const float support, const bool doUpres, 
                 const float &srcSize, const float &tgtSize);

  //--------------------------------------------------------------------------//

  V3f getDist(const V3i &doUpres, const V3f &srcP, const V3f &tgtP, 
              const V3f &srcSize, const V3f &tgtSize);

  //--------------------------------------------------------------------------//

  float getDist(const bool doUpres, const float &srcP, const float &tgtP, 
                const float &srcSize, const float &tgtSize);

  //--------------------------------------------------------------------------//

  template <typename Field_T, typename FilterOp_T, bool IsAnalytic_T>
  void separable(const Field_T &src, Field_T &tgt, const V3i &newRes,
                 const FilterOp_T &filterOp, const size_t dim)
  {
    typedef typename Field_T::value_type T;

    const V3i   srcRes    = src.dataWindow().size() + V3i(1);
    const float srcDomain = V3f(srcRes)[dim];
    const float tgtDomain = V3f(newRes)[dim];
    const float srcSize   = 1.0 / srcDomain;
    const float tgtSize   = 1.0 / tgtDomain;

    // Filter info
    const float support = filterOp.support();

    // Check if we're up-res'ing
    const bool doUpres = newRes[dim] > srcRes[dim] ? 1 : 0;

    // Resize the target
    tgt.setSize(newRes);

    // For each output voxel
    for (int k = 0; k < newRes.z; ++k) {
      for (int j = 0; j < newRes.y; ++j) {
        for (int i = 0; i < newRes.x; ++i) {
          T accumValue(filterOp.initialValue());
          if (IsAnalytic_T) {
            // Current position in target coordinates
            const float tgtP = discToCont(V3i(i, j ,k)[dim]);
            // Transform support to source coordinates
            std::pair<int, int> srcInterval = 
              srcSupportBBox(tgtP, support, doUpres, srcSize, tgtSize);
            // Clip against new data window
            srcInterval.first = 
              std::max(srcInterval.first, src.dataWindow().min[dim]);
            srcInterval.second = 
              std::min(srcInterval.second, src.dataWindow().max[dim]);
            // For each input voxel
            for (int s = srcInterval.first; s <= srcInterval.second; ++s) {
              // Index
              const int xIdx = dim == 0 ? s : i;
              const int yIdx = dim == 1 ? s : j;
              const int zIdx = dim == 2 ? s : k;
              // Value
              const T value = src.fastValue(xIdx, yIdx, zIdx);
              // Weights
              const float srcP   = discToCont(V3i(xIdx, yIdx, zIdx)[dim]);
              const float dist   = getDist(doUpres, srcP, tgtP, srcSize, tgtSize);
              const float weight = filterOp.eval(dist);
              // Update
              if (weight > 0.0f) {
                FilterOp_T::op(accumValue, value);
              }
            }
            // Update final value
            if (accumValue != static_cast<T>(filterOp.initialValue())) {
              tgt.fastLValue(i, j, k) = accumValue;
            }
          } else {
            float accumWeight = 0.0f;
            // Current position in target coordinates
            const float tgtP = discToCont(V3i(i, j ,k)[dim]);
            // Transform support to source coordinates
            std::pair<int, int> srcInterval = 
              srcSupportBBox(tgtP, support, doUpres, srcSize, tgtSize);
            // Clip against new data window
            srcInterval.first = 
              std::max(srcInterval.first, src.dataWindow().min[dim]);
            srcInterval.second = 
              std::min(srcInterval.second, src.dataWindow().max[dim]);
            // For each input voxel
            for (int s = srcInterval.first; s <= srcInterval.second; ++s) {
              // Index
              const int xIdx = dim == 0 ? s : i;
              const int yIdx = dim == 1 ? s : j;
              const int zIdx = dim == 2 ? s : k;
              // Value
              const T value      = src.fastValue(xIdx, yIdx, zIdx);
              // Weights
              const float srcP   = discToCont(V3i(xIdx, yIdx, zIdx)[dim]);
              const float dist   = getDist(doUpres, srcP, tgtP, srcSize, tgtSize);
              const float weight = filterOp.eval(dist);
              // Update
              accumWeight += weight;
              accumValue  += value * weight;
            }
            // Update final value
            if (accumWeight > 0.0f && accumValue != static_cast<T>(0.0)) {
              tgt.fastLValue(i, j, k) = accumValue / accumWeight;
            }
          }
        }
      }
    }
  }

  //--------------------------------------------------------------------------//

  //! Constant size for all dense fields
  template <typename Data_T>
  size_t threadingBlockSizeResample(const DenseField<Data_T> & /* f */)
  {
    return 16;
  }
  
  //! Use block size for sparse fields
  template <typename Data_T>
  size_t threadingBlockSizeResample(const SparseField<Data_T> &f)
  {
    return f.blockSize();
  }

  //--------------------------------------------------------------------------//

  template <typename Data_T>
  bool checkInputEmptyResample(const SparseField<Data_T> &src, 
                               const SparseField<Data_T> &/*tgt*/, 
                               const Box3i &tgtBox, const float support,
                                size_t dim)
  {
    const int intSupport = static_cast<int>(std::ceil(support * 0.5));
    const int pad        = std::max(0, intSupport);
    Box3i     tgtBoxPad  = tgtBox;
    tgtBoxPad.min[dim]  -= pad;
    tgtBoxPad.max[dim]  += pad;
    Box3i     srcBoxPad  = tgtBoxPad;
    // srcBoxPad.min[dim]  *= 2;
    // srcBoxPad.max[dim]  *= 2;

    // Get the block coordinates
    const Box3i dbsBounds = blockCoords(clipBounds(srcBoxPad, src.dataWindow()),
                                        &src);

    static boost::mutex mutex;
    boost::mutex::scoped_lock lock(mutex);

    // Check all blocks
    for (int k = dbsBounds.min.z; k <= dbsBounds.max.z; ++k) {
      for (int j = dbsBounds.min.y; j <= dbsBounds.max.y; ++j) {
        for (int i = dbsBounds.min.x; i <= dbsBounds.max.x; ++i) {
          if (src.blockIsAllocated(i, j, k) ||
              src.getBlockEmptyValue(i, j, k) != static_cast<Data_T>(0)) {
            return false;
          }
        }
      } 
    }

    // No hits. Empty
    return true;
  }

  //--------------------------------------------------------------------------//

  //! Fallback version always returns false
  template <typename Field_T>
  bool checkInputEmptyResample(const Field_T &/*src*/, const Field_T &/*tgt*/, 
                               const Box3i &/*tgtBox*/, const float /*support*/,
                               const size_t /*dim*/)
  {
    return false;
  }

  //--------------------------------------------------------------------------//

  template <typename Field_T, typename FilterOp_T, bool IsAnalytic_T>
  struct SeparableThreadOp
  {
    typedef typename Field_T::value_type T;

    SeparableThreadOp(const Field_T &src, Field_T &tgt, 
                      const FilterOp_T &filterOp, 
                      const size_t dim, 
                      const std::vector<Box3i> &blocks,
                      size_t &nextIdx, boost::mutex &mutex)
      : m_src(src),
        m_tgt(tgt),
        m_filterOp(filterOp), 
        m_dim(dim),
        m_blocks(blocks),
        m_nextIdx(nextIdx),
        m_mutex(mutex),
        m_numBlocks(blocks.size())
    {
    }

    void operator() () {
      typedef typename Field_T::value_type T;

      const V3i   srcRes    = m_src.dataWindow().size() + V3i(1);
      const V3i   newRes    = m_tgt.dataWindow().size() + V3i(1);
      const float srcDomain = V3f(srcRes)[m_dim];
      const float tgtDomain = V3f(newRes)[m_dim];
      const float srcSize   = 1.0 / srcDomain;
      const float tgtSize   = 1.0 / tgtDomain;

      // Filter info
      const float support = m_filterOp.support();

      // Check if we're up-res'ing
      const bool doUpres = newRes[m_dim] > srcRes[m_dim] ? 1 : 0;

      // Get next index to process
      size_t idx;
      {
        boost::mutex::scoped_lock lock(m_mutex);
        idx = m_nextIdx;
        m_nextIdx++;
      }
      // Keep going while there is data to process
      while (idx < m_numBlocks) {
        // Grab the bounds
        const Box3i box =  m_blocks[idx];
        // Early exit if input blocks are all empty
        if (!detail::checkInputEmptyResample(m_src, m_tgt, box, support, m_dim)) {
        // if (true) {
          // For each output voxel
          for (int k = box.min.z; k <= box.max.z; ++k) {
            for (int j = box.min.y; j <= box.max.y; ++j) {
              for (int i = box.min.x; i <= box.max.x; ++i) {
                T accumValue(m_filterOp.initialValue());
                if (IsAnalytic_T) {
                  // Current position in target coordinates
                  const float tgtP = discToCont(V3i(i, j ,k)[m_dim]);
                  // Transform support to source coordinates
                  std::pair<int, int> srcInterval = 
                    srcSupportBBox(tgtP, support, doUpres, srcSize, tgtSize);
                  // Clip against new data window
                  srcInterval.first = 
                    std::max(srcInterval.first, m_src.dataWindow().min[m_dim]);
                  srcInterval.second = 
                    std::min(srcInterval.second, m_src.dataWindow().max[m_dim]);
                  // For each input voxel
                  for (int s = srcInterval.first; s <= srcInterval.second; ++s) {
                    // Index
                    const int xIdx = m_dim == 0 ? s : i;
                    const int yIdx = m_dim == 1 ? s : j;
                    const int zIdx = m_dim == 2 ? s : k;
                    // Value
                    const T value = m_src.fastValue(xIdx, yIdx, zIdx);
                    // Weights
                    const float srcP   = discToCont(V3i(xIdx, yIdx, zIdx)[m_dim]);
                    const float dist   = getDist(doUpres, srcP, tgtP, srcSize, tgtSize);
                    const float weight = m_filterOp.eval(dist);
                    // Update
                    if (weight > 0.0f) {
                      FilterOp_T::op(accumValue, value);
                    }
                  }
                  // Update final value
                  if (accumValue != static_cast<T>(m_filterOp.initialValue())) {
                    m_tgt.fastLValue(i, j, k) = accumValue;
                  }
                } else {
                  float accumWeight = 0.0f;
                  // Current position in target coordinates
                  const float tgtP = discToCont(V3i(i, j ,k)[m_dim]);
                  // Transform support to source coordinates
                  std::pair<int, int> srcInterval = 
                    srcSupportBBox(tgtP, support, doUpres, srcSize, tgtSize);
                  // Clip against new data window
                  srcInterval.first = 
                    std::max(srcInterval.first, m_src.dataWindow().min[m_dim]);
                  srcInterval.second = 
                    std::min(srcInterval.second, m_src.dataWindow().max[m_dim]);
                  // For each input voxel
                  for (int s = srcInterval.first; s <= srcInterval.second; ++s) {
                    // Index
                    const int xIdx = m_dim == 0 ? s : i;
                    const int yIdx = m_dim == 1 ? s : j;
                    const int zIdx = m_dim == 2 ? s : k;
                    // Value
                    const T value      = m_src.fastValue(xIdx, yIdx, zIdx);
                    // Weights
                    const float srcP   = discToCont(V3i(xIdx, yIdx, zIdx)[m_dim]);
                    const float dist   = getDist(doUpres, srcP, tgtP, srcSize, tgtSize);
                    const float weight = m_filterOp.eval(dist);
                    // Update
                    accumWeight += weight;
                    accumValue  += value * weight;
                  }
                  // Update final value
                  if (accumWeight > 0.0f && accumValue != static_cast<T>(0.0)) {
                    m_tgt.fastLValue(i, j, k) = accumValue / accumWeight;
                  }
                }
              }
            }
          }
        }
        // Get next index
        {
          boost::mutex::scoped_lock lock(m_mutex);
          idx = m_nextIdx;
          m_nextIdx++;
        }
      }
    }

  private:

    // Data members ---

    const Field_T            &m_src;
    Field_T                  &m_tgt;
    const FilterOp_T         &m_filterOp;
    const size_t              m_dim;
    const std::vector<Box3i> &m_blocks;
    size_t                   &m_nextIdx;
    boost::mutex             &m_mutex;
    const size_t              m_numBlocks;
    
  };

  //--------------------------------------------------------------------------//

  template <typename Field_T, typename FilterOp_T, bool IsAnalytic_T>
  void separableThreaded(const Field_T &src, Field_T &tgt, const V3i &newRes,
                         const FilterOp_T &filterOp, const size_t dim,
                         const size_t numThreads)
  {
    // Resize the target
    tgt.setSize(newRes);

    // Determine granularity
    const size_t blockSize = threadingBlockSizeResample(src);

    // Build block list
    std::vector<Box3i> blocks;
    for (int k = 0; k < newRes.z; k += blockSize) {
      for (int j = 0; j < newRes.y; j += blockSize) {
        for (int i = 0; i < newRes.x; i += blockSize) {
          Box3i box;
          // Initialize block size
          box.min = V3i(i, j, k);
          box.max = box.min + V3i(blockSize - 1);
          // Clip against resolution
          box.max.x = std::min(box.max.x, newRes.x - 1);
          box.max.y = std::min(box.max.y, newRes.y - 1);
          box.max.z = std::min(box.max.z, newRes.z - 1);
          // Add to list
          blocks.push_back(box);
        }
      }
    }

    // Next index counter and mutex
    size_t nextIdx = 0;
    boost::mutex mutex;

    // Launch threads ---

    boost::thread_group threads;

    for (size_t i = 0; i < numThreads; ++i) {
      threads.create_thread(
        SeparableThreadOp<Field_T, FilterOp_T, FilterOp_T::isAnalytic >
        (src, tgt, filterOp, dim, blocks, nextIdx, mutex));
    }

    // Join
    threads.join_all();
  }

  //--------------------------------------------------------------------------//

  //! Resamples the source field into the target field, using separable
  //! execution, which is faster than resample().
  //! \note The extents of the field will be reset to match the data window.
  template <typename Field_T, typename FilterOp_T>
  bool separableResample(const Field_T &src, Field_T &tgt, const V3i &newRes,
                         const FilterOp_T &filterOp, const size_t numThreads)
  {
    using namespace detail;
  
    // typedef typename Field_T::value_type T;

    if (!src.dataWindow().hasVolume()) {
      return false;
    }

    if (src.dataWindow().min != V3i(0)) {
      return false;
    }

    // Temporary field for y component
    Field_T tmp;

    // Cache the old resolution
    V3i oldRes = src.dataWindow().size() + V3i(1);
    V3i xRes(newRes.x, oldRes.y, oldRes.z);
    V3i yRes(newRes.x, newRes.y, oldRes.z);
    V3i zRes(newRes.x, newRes.y, newRes.z);

    if (numThreads > 1) {
      // X axis (src into tgt)
      separableThreaded<Field_T, FilterOp_T, FilterOp_T::isAnalytic>(src, tgt, xRes, filterOp, 0, numThreads);
      // Y axis (tgt into temp)
      separableThreaded<Field_T, FilterOp_T, FilterOp_T::isAnalytic>(tgt, tmp, yRes, filterOp, 1, numThreads);
      // Z axis (temp into tgt)
      separableThreaded<Field_T, FilterOp_T, FilterOp_T::isAnalytic>(tmp, tgt, zRes, filterOp, 2, numThreads);
    } else {
      // X axis (src into tgt)
      separable<Field_T, FilterOp_T, FilterOp_T::isAnalytic>(src, tgt, xRes, filterOp, 0);
      // Y axis (tgt into temp)
      separable<Field_T, FilterOp_T, FilterOp_T::isAnalytic>(tgt, tmp, yRes, filterOp, 1);
      // Z axis (temp into tgt)
      separable<Field_T, FilterOp_T, FilterOp_T::isAnalytic>(tmp, tgt, zRes, filterOp, 2);
    }

    // Update final target with mapping and metadata
    tgt.name      = src.name;
    tgt.attribute = src.attribute;
    tgt.setMapping(src.mapping());
    tgt.copyMetadata(src);

    return true;
  }

  //--------------------------------------------------------------------------//

} // namespace detail

//----------------------------------------------------------------------------//
// Resizing function implementations
//----------------------------------------------------------------------------//

template <typename Field_T, typename FilterOp_T>
bool resample(const Field_T &src, Field_T &tgt, const V3i &newRes,
              const FilterOp_T &filterOp, const size_t numThreads)
{
  return detail::separableResample(src, tgt, newRes, filterOp, numThreads);
}

//----------------------------------------------------------------------------//

FIELD3D_NAMESPACE_HEADER_CLOSE

//----------------------------------------------------------------------------//

#endif // Include guard

