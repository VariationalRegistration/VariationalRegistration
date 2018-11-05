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

#include "itkVariationalRegistrationFilter.h"
#include "itkVariationalRegistrationMultiResolutionFilter.h"
#include "itkVariationalRegistrationDemonsFunction.h"
#include "itkVariationalRegistrationDiffusionRegularizer.h"
#include "itkVariationalRegistrationStopCriterion.h"
#include "itkVariationalRegistrationLogger.h"
#include "itkContinuousBorderWarpImageFilter.h"

#include "itkNearestNeighborInterpolateImageFunction.h"
#include "itkCommand.h"
#include "itkVectorCastImageFilter.h"
#include "itkImageFileWriter.h"


namespace{
// The following class is used to support callbacks
// on the filter in the pipeline that follows later
template<typename TRegistration>
class ShowProgressObject
{
public:
  ShowProgressObject(TRegistration* o)
    {m_Process = o;}
  void ShowProgress()
    {
    std::cout << "Progress: " << m_Process->GetProgress() << "  ";
    std::cout << "Iter: " << m_Process->GetElapsedIterations() << "  ";
    std::cout << "Metric: "   << m_Process->GetMetric()   << "  ";
    std::cout << "RMSChange: " << m_Process->GetRMSChange() << "  ";
    std::cout << std::endl;
    if ( m_Process->GetElapsedIterations() == 150 )
      { m_Process->StopRegistration(); }
    }
  typename TRegistration::Pointer m_Process;
};
}

// Template function to fill in an image with a circle.
template <typename TImage>
void
FillWithCircle(
TImage * image,
double * center,
double radius,
typename TImage::PixelType foregnd,
typename TImage::PixelType backgnd )
{

  using Iterator = itk::ImageRegionIteratorWithIndex<TImage>;
  Iterator it( image, image->GetBufferedRegion() );
  it.GoToBegin();

  typename TImage::IndexType index;
  double r2 = itk::Math::sqr( radius );

  while( !it.IsAtEnd() )
    {
    index = it.GetIndex();
    double distance = 0;
    for( unsigned int j = 0; j < TImage::ImageDimension; j++ )
      {
      distance += itk::Math::sqr((double) index[j] - center[j]);
      }
    if( distance <= r2 ) it.Set( foregnd );
    else it.Set( backgnd );
    ++it;
    }

}

// Template function to copy image regions
template <typename TImage>
void
CopyImageBuffer(
TImage *input,
TImage *output )
{
  using Iterator = itk::ImageRegionIteratorWithIndex<TImage>;
  Iterator outIt( output, output->GetBufferedRegion() );
  for( Iterator inIt( input, output->GetBufferedRegion() ); !inIt.IsAtEnd(); ++inIt, ++outIt )
    {
    outIt.Set( inIt.Get() );
    }

}

int VariationalRegistrationMultiResolutionFilterTest(int, char* [] )
{

  using PixelType = unsigned char;
  enum {ImageDimension = 2};
  using ImageType = itk::Image<PixelType,ImageDimension>;
  using VectorType = itk::Vector<float,ImageDimension>;
  using FieldType = itk::Image<VectorType,ImageDimension>;
  //using FloatImageType = itk::Image<VectorType::ValueType,ImageDimension>;
  using IndexType = ImageType::IndexType;
  using SizeType = ImageType::SizeType;
  using RegionType = ImageType::RegionType;

  //--------------------------------------------------------
  std::cout << "Generate input images and initial deformation field";
  std::cout << std::endl;

  ImageType::SizeValueType sizeArray[ImageDimension] = { 128, 128 };
  SizeType size;
  size.SetSize( sizeArray );

  IndexType index;
  index.Fill( 0 );

  RegionType region;
  region.SetSize( size );
  region.SetIndex( index );

  ImageType::Pointer moving = ImageType::New();
  ImageType::Pointer fixed = ImageType::New();
  FieldType::Pointer initField = FieldType::New();

  moving->SetLargestPossibleRegion( region );
  moving->SetBufferedRegion( region );
  moving->Allocate();

  fixed->SetLargestPossibleRegion( region );
  fixed->SetBufferedRegion( region );
  fixed->Allocate();

  initField->SetLargestPossibleRegion( region );
  initField->SetBufferedRegion( region );
  initField->Allocate();

  double center[ImageDimension];
  double radius;
  PixelType fgnd = 250;
  PixelType bgnd = 15;

  // fill moving with circle
  center[0] = 64; center[1] = 64; radius = 30;
  FillWithCircle<ImageType>( moving, center, radius, fgnd, bgnd );

  // fill fixed with circle
  center[0] = 62; center[1] = 64; radius = 32;
  FillWithCircle<ImageType>( fixed, center, radius, fgnd, bgnd );

  // fill initial deformation with zero vectors
  VectorType zeroVec;
  zeroVec.Fill( 0.0 );
  initField->FillBuffer( zeroVec );

  using CasterType = itk::VectorCastImageFilter<FieldType,FieldType>;
  CasterType::Pointer caster = CasterType::New();
  caster->SetInput( initField );
  caster->InPlaceOff();

  //-------------------------------------------------------------
  std::cout << "Run Multi-resolution registration and warp moving" << std::endl;

  // Setup registration function
  using DemonsFunctionType = itk::VariationalRegistrationDemonsFunction<
      ImageType, ImageType, FieldType>;
  DemonsFunctionType::Pointer demonsFunction = DemonsFunctionType::New();
  demonsFunction->SetGradientTypeToFixedImage();
  demonsFunction->SetTimeStep( 1.0 );

  // Setup regularizer
  using DiffusionRegularizerType = itk::VariationalRegistrationDiffusionRegularizer<FieldType>;
  DiffusionRegularizerType::Pointer diffRegularizer = DiffusionRegularizerType::New();
  diffRegularizer->SetAlpha( 0.1 );

  // Setup registration filter
  using RegistrationFilterType = itk::VariationalRegistrationFilter<
      ImageType,ImageType,FieldType>;
  RegistrationFilterType::Pointer regFilter = RegistrationFilterType::New();
  regFilter->SetRegularizer( diffRegularizer );
  regFilter->SetDifferenceFunction( demonsFunction );

  // Setup multi-resolution filter
  unsigned int its[2];
  its[1] = 200;
  its[0] = 200;

  using MRRegistrationFilterType = itk::VariationalRegistrationMultiResolutionFilter<ImageType,ImageType,FieldType>;
  MRRegistrationFilterType::Pointer mrRegFilter = MRRegistrationFilterType::New();
  mrRegFilter->SetRegistrationFilter( regFilter );
  mrRegFilter->SetMovingImage( moving );
  mrRegFilter->SetFixedImage( fixed );
  mrRegFilter->SetNumberOfLevels( 2 );
  mrRegFilter->SetNumberOfIterations( its );
  mrRegFilter->SetInitialField( caster->GetOutput() );

  // Setup stop criterion
  using StopCriterionType = itk::VariationalRegistrationStopCriterion<
      RegistrationFilterType,MRRegistrationFilterType>;
  StopCriterionType::Pointer stopCriterion = StopCriterionType::New();
  stopCriterion->SetRegressionLineSlopeThreshold( 0.0001 );
  stopCriterion->SetNumberOfFittingIterations( 20 );
  stopCriterion->SetMaxDistanceToRegressionLine(0.005);
  stopCriterion->LineFittingUseAbsoluteValuesOff();
  stopCriterion->PerformLineFittingMaxDistanceCheckOff();
  stopCriterion->PerformIncreaseCountCheckOn();
  stopCriterion->SetMultiResolutionPolicyToSimpleGraduated();

  regFilter->AddObserver( itk::IterationEvent(), stopCriterion );
  mrRegFilter->AddObserver( itk::IterationEvent(), stopCriterion );
  mrRegFilter->AddObserver( itk::InitializeEvent(), stopCriterion );

  // Setup logger
  using LoggerType = itk::VariationalRegistrationLogger<
      RegistrationFilterType,MRRegistrationFilterType>;
  LoggerType::Pointer logger = LoggerType::New();

  regFilter->AddObserver( itk::IterationEvent(), logger );
  mrRegFilter->AddObserver( itk::IterationEvent(), logger );

  // warp moving image
  using WarperType = itk::ContinuousBorderWarpImageFilter<ImageType,ImageType,FieldType>;
  WarperType::Pointer warper = WarperType::New();

  using CoordRepType = WarperType::CoordRepType;
  using InterpolatorType = itk::NearestNeighborInterpolateImageFunction<ImageType,CoordRepType>;
  InterpolatorType::Pointer interpolator = InterpolatorType::New();


  warper->SetInput( moving );
  warper->SetDisplacementField( mrRegFilter->GetOutput() );
  warper->SetInterpolator( interpolator );
  warper->SetOutputSpacing( fixed->GetSpacing() );
  warper->SetOutputOrigin( fixed->GetOrigin() );
  warper->SetOutputDirection( fixed->GetDirection() );
  warper->SetEdgePaddingValue( bgnd );

  warper->Print( std::cout );

  warper->Update();

  // ---------------------------------------------------------
  std::cout << "Compare warped moving and fixed." << std::endl;

  // compare the warp and fixed images
  itk::ImageRegionIterator<ImageType> fixedIter( fixed,
      fixed->GetBufferedRegion() );
  itk::ImageRegionIterator<ImageType> warpedIter( warper->GetOutput(),
      fixed->GetBufferedRegion() );

  unsigned int numPixelsDifferent = 0;
  while( !fixedIter.IsAtEnd() )
    {
    if( fixedIter.Get() != warpedIter.Get() )
      {
      numPixelsDifferent++;
      }
    ++fixedIter;
    ++warpedIter;
    }

  std::cout << "Number of pixels different: " << numPixelsDifferent;
  std::cout << std::endl;

  if( numPixelsDifferent > 10 )
    {
    std::cout << "Test failed - too many pixels different." << std::endl;
    return EXIT_FAILURE;
    }

  mrRegFilter->Print( std::cout );

  //--------------------------------------------------------------

  std::cout << "Test passed" << std::endl;
  return EXIT_SUCCESS;

}
