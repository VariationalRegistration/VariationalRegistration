/*=========================================================================
 *
 *  Copyright Insight Software Consortium
 *
 *  Licensed under the Apache License, Version 2.0 (the "License");
 *  you may not use this file except in compliance with the License.
 *  You may obtain a copy of the License at
 *
 *         http://www.apache.org/licenses/LICENSE-2.0.txt
 *
 *  Unless required by applicable law or agreed to in writing, software
 *  distributed under the License is distributed on an "AS IS" BASIS,
 *  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 *  See the License for the specific language governing permissions and
 *  limitations under the License.
 *
 *=========================================================================*/
#ifndef __itkVariationalRegistrationDiffusionRegularizer_h
#define __itkVariationalRegistrationDiffusionRegularizer_h

#include "itkVariationalRegistrationRegularizer.h"
#include "itkMultiThreader.h"

namespace itk {

/** \class itk::VariationalRegistrationDiffusionRegularizer
 *
 *  \sa VariationalRegistrationRegularizer
 *
 *  \ingroup VariationalRegistration
 */
template< class TDisplacementField >
class ITK_EXPORT VariationalRegistrationDiffusionRegularizer
  : public VariationalRegistrationRegularizer< TDisplacementField >
{
public:
  /** Standard class typedefs */
  typedef VariationalRegistrationDiffusionRegularizer  Self;
  typedef VariationalRegistrationRegularizer<
      TDisplacementField >                       Superclass;
  typedef SmartPointer< Self >                   Pointer;
  typedef SmartPointer< const Self >             ConstPointer;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** Run-time type information (and related methods) */
  itkTypeMacro(itkVariationalRegistrationDiffusionRegularizer, itkVariationalRegistrationRegularizer);

  /** Dimensionality of input and output data is assumed to be the same. */
  itkStaticConstMacro(ImageDimension, unsigned int, TDisplacementField::ImageDimension);

  /** Deformation field types, inherited from Superclass. */
  typedef typename Superclass::DisplacementFieldType         DisplacementFieldType;
  typedef typename Superclass::DisplacementFieldPointer      DisplacementFieldPointer;
  typedef typename Superclass::DisplacementFieldConstPointer DisplacementFieldConstPointer;
  typedef typename Superclass::PixelType                     PixelType;

  typedef typename Superclass::ValueType                     ValueType;

  /** Types for buffer image. */
  typedef Image<ValueType, ImageDimension>                   BufferImageType;
  typedef typename BufferImageType::Pointer                  BufferImagePointer;
  typedef typename BufferImageType::RegionType               BufferImageRegionType;

  /** Set/Get the regularization weight alpha */
  itkSetMacro( Alpha, float );
  itkGetMacro( Alpha, float );

protected:
  VariationalRegistrationDiffusionRegularizer();
  ~VariationalRegistrationDiffusionRegularizer() {}
  void PrintSelf(std::ostream& os, Indent indent) const;

  /** Execute regularization. This method is multi-threaded but does not
   * use ThreadedGenerateData(). */
  virtual void GenerateData();

  /** Method for initialization. Buffer images are allocated and the matrices
   * calculated in this method. */
  virtual void Initialize();

  /** Calculation and LU decomposition of the tridiagonal matrices
   * using the Thomas algorithm (TDMA). */
  void InitLUMatrices( ValueType** alpha, ValueType** beta, ValueType** gamma,
      int n, int dim );

  /** Regularize a given component (dimension) of the deformation field. This
   *  is called for each ImageDimension by GenerateData(). */
  void RegularizeComponent( const int component );

  /** A struct to store parameters for multithreaded function call. */
  struct CalcBufferThreadStruct
  {
    VariationalRegistrationDiffusionRegularizer *Filter;
    unsigned int component;        // The current dimension.
    BufferImagePointer bPtr;       // Pointer to the image buffer.
  };

  /** A struct to store parameters for multithreaded function call. */
  struct RegularizeThreadStruct
  {
    VariationalRegistrationDiffusionRegularizer *Filter;
    unsigned int direction;        // Current direction.
    ValueType* alpha;              // Pointer to matrix diagonal.
    ValueType* beta;               // Pointer to matrix subdiagonal.
    ValueType* gamma;              // Pointer to matrix superdiagonal.
    BufferImagePointer bPtr;       // Pointer to force field image buffer.
    BufferImagePointer vPtr;       // Pointer to temporal result image buffer.
  };

  /** A struct to store parameters for multithreaded function call. */
  struct MergeDirectionsThreadStruct
  {
    VariationalRegistrationDiffusionRegularizer *Filter;
    unsigned int component;        // The current dimension.
    BufferImagePointer* vPtr;      // Pointer to temporal image buffers.
  };

  /** Method for multi-threaded calculation of the image buffer. */
  static ITK_THREAD_RETURN_TYPE CalcBufferCallback( void *arg );

  /** Method for multi-threaded regularization of the image buffer. */
  static ITK_THREAD_RETURN_TYPE RegularizeDirectionCallback( void *arg );

  /** Method for multi-threaded calculation of the final field. */
  static ITK_THREAD_RETURN_TYPE MergeDirectionsCallback( void *arg );

  /** Split the boundary face orthogonal to "inDir" into "num" pieces, returning
   * region "i" as "splitRegion". This method is called "num" times. The
   * regions must not overlap. The method returns the number of pieces that
   * the routine is capable of splitting the output RequestedRegion,
   * i.e. return value is less than or equal to "num". */
  virtual int SplitBoundaryFaceRegion( int i, int num, int inDir,
      BufferImageRegionType& splitRegion );

private:
  VariationalRegistrationDiffusionRegularizer(const Self&); //purposely not implemented
  void operator=(const Self&); //purposely not implemented

  /** Weight of the regularization term. */
  float m_Alpha;

  /** The size of the displacement field. */
  typename DisplacementFieldType::SizeType m_Size;

  /** The spacing of the displacement field. */
  typename DisplacementFieldType::SpacingType m_Spacing;

  // Attributes for AOS calculation
  /** Pointer to a temporal image for the regularization.  */
  BufferImagePointer m_BufferImage;

  /** Buffers for the regularized fields in each direction. Stored in separate
   * images instead of a vector image for better memory management. */
  BufferImagePointer m_V[ImageDimension];

  /** Array for the diagonals of the factorized matrices for each dimension */
  ValueType* m_MatrixAlpha[ImageDimension];

  /** Array for the subdiagonals of the factorized matrices for each dimension */
  ValueType* m_MatrixBeta[ImageDimension];

  /** Array for the superdiagonals of the factorized matrices for each dimension */
  ValueType* m_MatrixGamma[ImageDimension];
};

}

#ifndef ITK_MANUAL_INSTANTIATION
# include "itkVariationalRegistrationDiffusionRegularizer.hxx"
#endif

#endif