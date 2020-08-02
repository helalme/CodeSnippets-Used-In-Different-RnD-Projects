/*
 * Test filter for generating lamellae along 10 different directions
*/

#include "GenerateLamellaeForDuplex.h"

#include "SIMPLib/Common/Constants.h"

#include "SIMPLib/FilterParameters/AbstractFilterParametersReader.h"
#include "SIMPLib/FilterParameters/DataArrayCreationFilterParameter.h"
#include "SIMPLib/FilterParameters/DataArraySelectionFilterParameter.h"
#include "SIMPLib/FilterParameters/LinkedChoicesFilterParameter.h"
#include "SIMPLib/FilterParameters/SeparatorFilterParameter.h"
#include "SIMPLib/FilterParameters/StringFilterParameter.h"
#include "SIMPLib/Geometry/ImageGeom.h"

#include "Generic/GenericConstants.h"
#include "Generic/GenericVersion.h"
#include "OrientationAnalysis/OrientationAnalysisConstants.h"
#include "OrientationAnalysis/OrientationAnalysisVersion.h"
#include <stdlib.h>

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
GenerateLamellaeForDuplex::GenerateLamellaeForDuplex()
: m_FeatureIdsArrayPath(SIMPL::Defaults::ImageDataContainerName, SIMPL::Defaults::CellAttributeMatrixName, SIMPL::CellData::FeatureIds)
, m_CellPhasesArrayPath(SIMPL::Defaults::ImageDataContainerName, SIMPL::Defaults::CellAttributeMatrixName, SIMPL::CellData::Phases)
//, m_FeaturePhasesArrayPath(SIMPL::Defaults::ImageDataContainerName, SIMPL::Defaults::CellFeatureAttributeMatrixName, SIMPL::FeatureData::Phases)
, m_FeatureRectArrayPath(SIMPL::Defaults::ImageDataContainerName, SIMPL::Defaults::CellFeatureAttributeMatrixName, SIMPL::FeatureData::FeatureRect)
, m_EulerAnglesArrayPath(SIMPL::Defaults::ImageDataContainerName, SIMPL::Defaults::CellAttributeMatrixName, SIMPL::CellData::EulerAngles)
, m_AvgEulerAnglesArrayPath(SIMPL::Defaults::ImageDataContainerName, SIMPL::Defaults::CellFeatureAttributeMatrixName, SIMPL::FeatureData::AvgEulerAngles)
, m_FeatureIds(nullptr)
, m_CellPhases(nullptr)
//, m_FeaturePhases(nullptr)
, m_FeatureRect(nullptr)
, m_EulerAngles(nullptr)
, m_AvgEulerAngles(nullptr)
{
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
GenerateLamellaeForDuplex::~GenerateLamellaeForDuplex() = default;

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void GenerateLamellaeForDuplex::initialize()
{
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void GenerateLamellaeForDuplex::setupFilterParameters()
{
  FilterParameterVector parameters;
  parameters.push_back(SeparatorFilterParameter::New("Cell Data", FilterParameter::RequiredArray));
  {
    DataArraySelectionFilterParameter::RequirementType req = DataArraySelectionFilterParameter::CreateRequirement(SIMPL::TypeNames::Int32, 1, AttributeMatrix::Type::Cell, IGeometry::Type::Image);
    parameters.push_back(SIMPL_NEW_DA_SELECTION_FP("Feature Ids", FeatureIdsArrayPath, FilterParameter::RequiredArray, GenerateLamellaeForDuplex, req));
  }
  {
    DataArraySelectionFilterParameter::RequirementType req = DataArraySelectionFilterParameter::CreateRequirement(SIMPL::TypeNames::Int32, 1, AttributeMatrix::Type::Cell, IGeometry::Type::Image);
    parameters.push_back(SIMPL_NEW_DA_SELECTION_FP("Phases", CellPhasesArrayPath, FilterParameter::RequiredArray, GenerateLamellaeForDuplex, req));
  }
  {
    DataArraySelectionFilterParameter::RequirementType req = DataArraySelectionFilterParameter::CreateRequirement(SIMPL::TypeNames::Float, 3, AttributeMatrix::Type::Cell, IGeometry::Type::Image);
    parameters.push_back(SIMPL_NEW_DA_SELECTION_FP("Euler Angles", EulerAnglesArrayPath, FilterParameter::RequiredArray, GenerateLamellaeForDuplex, req));
  }

  parameters.push_back(SeparatorFilterParameter::New("Cell Feature Data", FilterParameter::RequiredArray));
  {
    DataArraySelectionFilterParameter::RequirementType req =
        DataArraySelectionFilterParameter::CreateRequirement(SIMPL::TypeNames::Float, 3, AttributeMatrix::Type::CellFeature, IGeometry::Type::Image);
    parameters.push_back(SIMPL_NEW_DA_SELECTION_FP("Average Euler Angles", AvgEulerAnglesArrayPath, FilterParameter::RequiredArray, GenerateLamellaeForDuplex, req, 0));
  }

  parameters.push_back(SeparatorFilterParameter::New("Cell Feature Data", FilterParameter::CreatedArray));

  DataArrayCreationFilterParameter::RequirementType dacReq;
  dacReq.amTypes = {AttributeMatrix::Type::CellFeature};
  parameters.push_back(SIMPL_NEW_DA_CREATION_FP("Feature Rect", FeatureRectArrayPath, FilterParameter::CreatedArray, GenerateLamellaeForDuplex, dacReq));

  setFilterParameters(parameters);
}

void GenerateLamellaeForDuplex::readFilterParameters(AbstractFilterParametersReader* reader, int index)
{
  reader->openFilterGroup(this, index);
  // setFeaturePhasesArrayPath(reader->readDataArrayPath("FeaturePhasesArrayPath", getFeaturePhasesArrayPath()));
  setCellPhasesArrayPath(reader->readDataArrayPath("CellPhasesArrayPath", getCellPhasesArrayPath()));
  setFeatureIdsArrayPath(reader->readDataArrayPath("FeatureIdsArrayPath", getFeatureIdsArrayPath()));
  setAvgEulerAnglesArrayPath(reader->readDataArrayPath("AvgEulerAnglesArrayPath", getAvgEulerAnglesArrayPath()));
  setEulerAnglesArrayPath(reader->readDataArrayPath("EulerAnglesArrayPath", getEulerAnglesArrayPath()));
  setFeatureRectArrayPath(reader->readDataArrayPath("FeatureRectArrayPath", getFeatureRectArrayPath()));
  reader->closeFilterGroup();
}
// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void GenerateLamellaeForDuplex::dataCheck()
{
  setErrorCondition(0);
  setWarningCondition(0);

  QVector<size_t> cDims(1, 1);
  m_FeatureIdsPtr = getDataContainerArray()->getPrereqArrayFromPath<DataArray<int32_t>, AbstractFilter>(this, getFeatureIdsArrayPath(),
                                                                                                        cDims); /* Assigns the shared_ptr<> to an instance variable that is a weak_ptr<> */
  if(nullptr != m_FeatureIdsPtr.lock())                                                                         /* Validate the Weak Pointer wraps a non-nullptr pointer to a DataArray<T> object */
  {
    m_FeatureIds = m_FeatureIdsPtr.lock()->getPointer(0);
  }

  m_CellPhasesPtr = getDataContainerArray()->getPrereqArrayFromPath<DataArray<int32_t>, AbstractFilter>(this, getCellPhasesArrayPath(),
                                                                                                        cDims); /* Assigns the shared_ptr<> to an instance variable that is a weak_ptr<> */
  if(nullptr != m_CellPhasesPtr.lock())                                                                         /* Validate the Weak Pointer wraps a non-nullptr pointer to a DataArray<T> object */
  {
    m_CellPhases = m_CellPhasesPtr.lock()->getPointer(0);
  } /* Now assign the raw pointer to data from the DataArray<T> object */

  cDims[0] = 6;
  {
    m_FeatureRectPtr = getDataContainerArray()->createNonPrereqArrayFromPath<DataArray<size_t>, AbstractFilter, size_t>(
        this, getFeatureRectArrayPath(), 0, cDims); // Assigns the shared_ptr<> to an instance variable that is a weak_ptr<> /
    if(nullptr != m_FeatureRectPtr.lock())          // Validate the Weak Pointer wraps a non-nullptr pointer to a DataArray<T> object /
    {
      m_FeatureRect = m_FeatureRectPtr.lock()->getPointer(0);
    } // Now assign the raw pointer to data from the DataArray<T> object //
  }

  cDims[0] = 3;
  m_EulerAnglesPtr = getDataContainerArray()->getPrereqArrayFromPath<DataArray<float>, AbstractFilter>(this, getEulerAnglesArrayPath(),
                                                                                                       cDims); /* Assigns the shared_ptr<> to an instance variable that is a weak_ptr<> */
  if(nullptr != m_EulerAnglesPtr.lock())                                                                       /* Validate the Weak Pointer wraps a non-nullptr pointer to a DataArray<T> object */
  {
    m_EulerAngles = m_EulerAnglesPtr.lock()->getPointer(0);
  } /* Now assign the raw pointer to data from the DataArray<T> object */

  cDims[0] = 3;
  m_AvgEulerAnglesPtr = getDataContainerArray()->getPrereqArrayFromPath<DataArray<float>, AbstractFilter>(this, getAvgEulerAnglesArrayPath(),
                                                                                                          cDims); /* Assigns the shared_ptr<> to an instance variable that is a weak_ptr<> */
  if(nullptr != m_AvgEulerAnglesPtr.lock())                                                                       /* Validate the Weak Pointer wraps a non-nullptr pointer to a DataArray<T> object */
  {
    m_AvgEulerAngles = m_AvgEulerAnglesPtr.lock()->getPointer(0);
  } /* Now assign the raw pointer to data from the DataArray<T> object */
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void GenerateLamellaeForDuplex::preflight()
{
  // These are the REQUIRED lines of CODE to make sure the filter behaves correctly
  setInPreflight(true);              // Set the fact that we are preflighting.
  emit preflightAboutToExecute();    // Emit this signal so that other widgets can do one file update
  emit updateFilterParameters(this); // Emit this signal to have the widgets push their values down to the filter
  dataCheck();                       // Run our DataCheck to make sure everthing is setup correctly
  emit preflightExecuted();          // We are done preflighting this filter
  setInPreflight(false);             // Inform the system this filter is NOT in preflight mode anymore.
}
// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void GenerateLamellaeForDuplex::generate_lamellaeX()
{
  //*****************************finding feature corner***** begin*****************************
  int featureId = 0;
  size_t numComps = 6;
  QVector<size_t> cDims(1, numComps);
  int err = 0;
  AttributeMatrix::Pointer featureAM = getDataContainerArray()->getPrereqAttributeMatrixFromPath<AbstractFilter>(this, m_FeatureRectArrayPath, err);

  Int32ArrayType::Pointer cellFeatureIds = m_FeatureIdsPtr.lock();

  // Create corners array, which stores pixel coordinates for the top-left and bottom-right coordinates of each feature object
  UInt64ArrayType::Pointer corners = m_FeatureRectPtr.lock();
  for(size_t i = 0; i < corners->getNumberOfTuples(); i++)
  {
    corners->setComponent(i, 0, std::numeric_limits<uint32_t>::max());
    corners->setComponent(i, 1, std::numeric_limits<uint32_t>::max());
    corners->setComponent(i, 2, std::numeric_limits<uint32_t>::max());
    corners->setComponent(i, 3, std::numeric_limits<uint32_t>::min());
    corners->setComponent(i, 4, std::numeric_limits<uint32_t>::min());
    corners->setComponent(i, 5, std::numeric_limits<uint32_t>::min());
  }
  AttributeMatrix::Pointer featureIdsAM = getDataContainerArray()->getAttributeMatrix(m_FeatureIdsArrayPath);
  QVector<size_t> imageDims = featureIdsAM->getTupleDimensions();
  size_t xDim = imageDims[0], yDim = imageDims[1], zDim = imageDims[2];

  size_t index = 0;
  // Store the coordinates in the corners array
  for(uint32_t z = 0; z < zDim; z++)
  {
    if(getCancel())
    {
      return;
    }

    for(uint32_t y = 0; y < yDim; y++)
    {
      for(uint32_t x = 0; x < xDim; x++)
      {
        index = sub2ind(imageDims, x, y, z); // Index into cellFeatureIds array

        featureId = m_FeatureIds[index];
        if(featureId == 0)
        {
          continue;
        }

        if(featureId >= corners->getNumberOfTuples())
        {
          setErrorCondition(-31000);
          QString ss = QObject::tr("The feature attribute matrix '%1' has a smaller tuple count than the maximum feature id in '%2'").arg(featureAM->getName()).arg(cellFeatureIds->getName());
          notifyErrorMessage(getHumanLabel(), ss, getErrorCondition());
          return;
        }

        uint64_t* featureCorner = corners->getPointer(featureId * numComps);

        uint32_t val = featureCorner[0];
        if(x < featureCorner[0])
        {
          featureCorner[0] = x;
        }
        val = featureCorner[1];
        if(y < featureCorner[1])
        {
          featureCorner[1] = y;
        }
        val = featureCorner[2];
        if(z < featureCorner[2])
        {
          featureCorner[2] = z;
        }

        val = featureCorner[3];
        if(x > featureCorner[3])
        {
          featureCorner[3] = x;
        }
        val = featureCorner[4];
        if(y > featureCorner[4])
        {
          featureCorner[4] = y;
        }
        val = featureCorner[5];
        if(z > featureCorner[5])
        {
          featureCorner[5] = z;
        }
      }
    }
  }
  //*****************************finding feature corner***** end*****************************
  // Feature corner charao X, Y, Z directional lamellae generate kora jabe.
  // 45 degree lamellae jevabe generate kora hoise thik shevabe algorithm ready korte hobe
  // even 45 er ordhek direction eo lamellae generate kora jabe
  //*****************************Generating Lamellae***** start*****************************

  DataContainer::Pointer m = getDataContainerArray()->getDataContainer(getFeatureIdsArrayPath().getDataContainerName());
  ImageGeom::Pointer imageGeom = m->getGeometryAs<ImageGeom>();

  // Get number of voxels from imageGeometry along X, Y, Z directions
  size_t xPoints = imageGeom->getXPoints();
  size_t yPoints = imageGeom->getYPoints();
  size_t zPoints = imageGeom->getZPoints();

  // yStride and zStride are variables to track serial number of each voxel from xPoints, yPoints and zPoints
  size_t zStride = 0;
  size_t yStride = 0;
  // FeaturesWithXLamellae is to list features with X-plane lamellae, similarly FeaturesWithYLamellae and FeaturesWithZLamellae
  QVector<int32_t> FeaturesWithXLamellae;
  QVector<int32_t> FeaturesWithYLamellae;
  QVector<int32_t> FeaturesWithZLamellae;
  QVector<int32_t> FeaturesWithXY45degreeLamellae;
  QVector<int32_t> FeaturesWithXYminus45degreeLamellae;
  QVector<int32_t> FeaturesWithYZ45degreeLamellae;
  QVector<int32_t> FeaturesWithYZminus45degreeLamellae;
  QVector<int32_t> FeaturesWithZX45degreeLamellae;
  QVector<int32_t> FeaturesWithZXminus45degreeLamellae;
  QVector<int32_t> FeaturesWithXY60degreeLamellae;
  QVector<int32_t> RefCellForXY60degreeLamellae;

  /* the number of grains went for lamellae generation, e.g. count=1 means one grain is set for generating X-lamellae,
  count=2 means one grain is set for generating X-lamellae, second one is for Y lamellae,
  count=3 means one grain is set for generating X-lamellae, one for Y lamellae and one for Z lamellae
  */
  int32_t count = 0;
  // int32_t countXY = 0;

  // i, j, k for iterating voxel coordinates along z, y and x axes respectively
  for(size_t i = 0; i < zPoints; i++)
  {
    zStride = i * xPoints * yPoints;
    for(size_t j = 0; j < yPoints; j++)
    {
      yStride = j * xPoints;
      for(size_t k = 0; k < xPoints; k++)
      {
        // cellNumber is the serial number of a voxel. gnum is the FeatureId of that voxel (cellNumber)
        size_t cellNumber = zStride + yStride + k;
        int32_t gnum = m_FeatureIds[cellNumber];

        // QVector::contains(gnum) is similar to std::find(vector.begin(), vector.end(), gnum)
        /*main algorithm starts here: based on Voxel coordinate iterate each voxel till the end, check the featurePhase containing that
    voxel wheather it should be considered for lamellae generation or not. If yes, check whether that feature already considered for
    X or Y or Z lamellae generation. If yes, simply direct to the respective X or Y or Z block of if-elseif-else statement, if not
    go to the else part of if-elseif-else and specify which directional lamella the feature should have. Continue for the next voxel
    and complete till the last voxel  */
        if((int)m_CellPhases[cellNumber] == 1)     // only cellPhase=1 is being considered for lamellae generation, this is for Duplex type. 
        {							//fully lamellar korte hole ei if{} block ta commented kore dao
          if(FeaturesWithXLamellae.contains(gnum)) // is gnum already considered for X-Plane-lamellae?
          {
            size_t distanceFromXmin = k - corners->getComponent(gnum, 0); // how far the current voxels x-coordinate from min x of the current feature
            size_t lamellarNumber = distanceFromXmin / 5;                 // lamellar width is 5 here. koto tomo lamella in the current grain
            size_t phaseSelection = lamellarNumber % 2;                   // phase of the current lamella thik korte hobe  phaseSelection theke

            if(phaseSelection == 0)
            {
              m_CellPhases[cellNumber] = 1;
              m_EulerAngles[cellNumber * 3] = m_EulerAngles[gnum * 3] + 1.0; // Feature EulerAngle ke Cell e transfer kora hosse 1 jog kore
            }
            else
            {
              m_CellPhases[cellNumber] = 2;
              m_EulerAngles[cellNumber * 3] = m_EulerAngles[gnum * 3] - 1.0;
              if(m_EulerAngles[cellNumber * 3] < 0.0)
                m_EulerAngles[cellNumber * 3] = 0.0;
            }
          }
          else if(FeaturesWithYLamellae.contains(gnum)) // is gnum already considered for Y-Plane-lamellae?
          {
            size_t distanceFromYmin = j - corners->getComponent(gnum, 1);
            size_t lamellarNumber = distanceFromYmin / 5;
            size_t phaseSelection = lamellarNumber % 2;

            if(phaseSelection == 0)
            {
              m_CellPhases[cellNumber] = 1;
              m_EulerAngles[cellNumber * 3 + 1] = m_EulerAngles[gnum * 3 + 1] + 1.0;
            }
            else
            {
              m_CellPhases[cellNumber] = 2;
              m_EulerAngles[cellNumber * 3 + 1] = m_EulerAngles[gnum * 3 + 1] - 1.0;
              if(m_EulerAngles[cellNumber * 3 + 1] < 0.0)
                m_EulerAngles[cellNumber * 3 + 1] = 0.0;
            }
          }
          else if(FeaturesWithZLamellae.contains(gnum)) // is gnum already considered for Z-Plane-lamellae?
          {
            size_t distanceFromZmin = i - corners->getComponent(gnum, 2);
            size_t lamellarNumber = distanceFromZmin / 5;
            size_t phaseSelection = lamellarNumber % 2;

            if(phaseSelection == 0)
            {
              m_CellPhases[cellNumber] = 1;
              m_EulerAngles[cellNumber * 3 + 2] = m_EulerAngles[gnum * 3 + 2] + 1.0;
            }
            else
            {
              m_CellPhases[cellNumber] = 2;
              m_EulerAngles[cellNumber * 3 + 2] = m_EulerAngles[gnum * 3 + 2] - 1.0;
              if(m_EulerAngles[cellNumber * 3 + 2] < 0.0)
                m_EulerAngles[cellNumber * 3 + 2] = 0.0;
            }
          }
          else if(FeaturesWithXY45degreeLamellae.contains(gnum)) // is gnum already considered for lamellae along XY Plane 45degree?
          {
            int refXY = (int)k - (int)j;
            if(lamellarDecisionXYplane(refXY) == true)
            {
              m_CellPhases[cellNumber] = 2;
              m_EulerAngles[cellNumber * 3] = m_EulerAngles[gnum * 3] + 1.0;
            }
            else
            {
              m_CellPhases[cellNumber] = 1;
              m_EulerAngles[cellNumber * 3] = m_EulerAngles[gnum * 3] - 1.0;
              if(m_EulerAngles[cellNumber * 3] < 0.0)
                m_EulerAngles[cellNumber * 3] = 0.0;
            }
          }
          else if(FeaturesWithXYminus45degreeLamellae.contains(gnum)) // is gnum already considered for lamellae along XY Plane minus45degree?
          {
            int refXY = (int)k + (int)j;
            if(lamellarDecisionXYplane(refXY) == true)
            {
              m_CellPhases[cellNumber] = 2;
              m_EulerAngles[cellNumber * 3] = m_EulerAngles[gnum * 3] + 1.0;
            }
            else
            {
              m_CellPhases[cellNumber] = 1;
              m_EulerAngles[cellNumber * 3] = m_EulerAngles[gnum * 3] - 1.0;
              if(m_EulerAngles[cellNumber * 3] < 0.0)
                m_EulerAngles[cellNumber * 3] = 0.0;
            }
          }
          else if(FeaturesWithYZ45degreeLamellae.contains(gnum)) // is gnum already considered for lamellae along YZ Plane 45degree?
          {
            int refYZ = (int)j - (int)i;
            if(lamellarDecisionXYplane(refYZ) == true)
            {
              m_CellPhases[cellNumber] = 2;
              m_EulerAngles[cellNumber * 3 + 1] = m_EulerAngles[gnum * 3 + 1] + 1.0;
            }
            else
            {
              m_CellPhases[cellNumber] = 1;
              m_EulerAngles[cellNumber * 3 + 1] = m_EulerAngles[gnum * 3 + 1] - 1.0;
              if(m_EulerAngles[cellNumber * 3 + 1] < 0.0)
                m_EulerAngles[cellNumber * 3 + 1] = 0.0;
            }
          }
          else if(FeaturesWithYZminus45degreeLamellae.contains(gnum)) // is gnum already considered for lamellae along YZ Plane minus45degree?
          {
            int refYZ = (int)j + (int)i;
            if(lamellarDecisionXYplane(refYZ) == true)
            {
              m_CellPhases[cellNumber] = 2;
              m_EulerAngles[cellNumber * 3 + 1] = m_EulerAngles[gnum * 3 + 1] + 1.0;
            }
            else
            {
              m_CellPhases[cellNumber] = 1;
              m_EulerAngles[cellNumber * 3 + 1] = m_EulerAngles[gnum * 3 + 1] - 1.0;
              if(m_EulerAngles[cellNumber * 3 + 1] < 0.0)
                m_EulerAngles[cellNumber * 3 + 1] = 0.0;
            }
          }
          else if(FeaturesWithZX45degreeLamellae.contains(gnum)) // is gnum already considered for lamellae along ZX Plane 45degree?
          {
            int refYZ = (int)i - (int)k;
            if(lamellarDecisionXYplane(refYZ) == true)
            {
              m_CellPhases[cellNumber] = 2;
              m_EulerAngles[cellNumber * 3 + 2] = m_EulerAngles[gnum * 3 + 2] + 1.0;
            }
            else
            {
              m_CellPhases[cellNumber] = 1;
              m_EulerAngles[cellNumber * 3 + 2] = m_EulerAngles[gnum * 3 + 2] - 1.0;
              if(m_EulerAngles[cellNumber * 3 + 2] < 0.0)
                m_EulerAngles[cellNumber * 3 + 2] = 0.0;
            }
          }
          else if(FeaturesWithZXminus45degreeLamellae.contains(gnum)) // is gnum already considered for lamellae along ZX Plane minus45degree?
          {
            int refYZ = (int)i + (int)k;
            if(lamellarDecisionXYplane(refYZ) == true)
            {
              m_CellPhases[cellNumber] = 2;
              m_EulerAngles[cellNumber * 3 + 2] = m_EulerAngles[gnum * 3 + 2] + 1.0;
            }
            else
            {
              m_CellPhases[cellNumber] = 1;
              m_EulerAngles[cellNumber * 3 + 2] = m_EulerAngles[gnum * 3 + 2] - 1.0;
              if(m_EulerAngles[cellNumber * 3 + 2] < 0.0)
                m_EulerAngles[cellNumber * 3 + 2] = 0.0;
            }
          }
          /*
          else if(FeaturesWithXY60degreeLamellae.contains(gnum))
          {
            auto it = FeaturesWithXY60degreeLamellae.indexOf(gnum);
            int refX = RefCellForXY60degreeLamellae[it * 2];
            int refY = RefCellForXY60degreeLamellae[it * 2 + 1];
            int delY = abs(int(j - refY));
            int levelRefX = refX + (delY + 1) / 2;
            int levelRefY = refY + delY;
            int delXfromLevelRefX = k - levelRefX;
            if(delXfromLevelRefX >= 0)
            {
            if((delXfromLevelRefX / 5) % 2 == 0)
            {
              m_CellPhases[cellNumber] = 1;
              m_EulerAngles[cellNumber * 3] = m_EulerAngles[gnum * 3] + 1.0;
            }
            else
            {
              m_CellPhases[cellNumber] = 2;
              m_EulerAngles[cellNumber * 3] = m_EulerAngles[gnum * 3] - 1.0;
              if(m_EulerAngles[cellNumber * 3] < 0.0)
              m_EulerAngles[cellNumber * 3] = 0.0;
            }
            }
          

            else
            {
            if(((abs(int(delXfromLevelRefX)) - 1) / 5) % 2 == 1)
            {
              m_CellPhases[cellNumber] = 1;
              m_EulerAngles[cellNumber * 3] = m_EulerAngles[gnum * 3] + 1.0;
            }
            else
            {
              m_CellPhases[cellNumber] = 2;
              m_EulerAngles[cellNumber * 3] = m_EulerAngles[gnum * 3] - 1.0;
              if(m_EulerAngles[cellNumber * 3] < 0.0)
              m_EulerAngles[cellNumber * 3] = 0.0;
            }
            }
          }
          */
          else // enter here if the voxel is not considered yet wheather it should be gone for lamellae or not. If yes, which direction?
          {
            count++;
            if(count % 10 == 1) // 1st, 4th, 7th.... grains with phase 2 will go for X-Plane-Lamellae, if lamellar directions are 3
            {
              FeaturesWithXLamellae.push_back(gnum);
              size_t distanceFromXmin = k - corners->getComponent(gnum, 0);
              size_t lamellarNumber = distanceFromXmin / 5;
              size_t phaseSelection = lamellarNumber % 2;
              // uint64_t midPoint = corners->getComponent(gnum, 0) + (corners->getComponent(gnum, 3) - corners->getComponent(gnum, 0)) / 2;

              if(phaseSelection == 0)
              {
                m_CellPhases[cellNumber] = 1;
                m_EulerAngles[cellNumber * 3] = m_EulerAngles[gnum * 3] + 1.0;
              }
              else
              {
                m_CellPhases[cellNumber] = 2;
                m_EulerAngles[cellNumber * 3] = m_EulerAngles[gnum * 3] - 1.0;
                if(m_EulerAngles[cellNumber * 3] < 0.0)
                  m_EulerAngles[cellNumber * 3] = 0.0;
              }
            }
            else if(count % 10 == 2) // 2nd, 5th, 8th.... grains with phase 2 will go for Y-Plane-Lamellae
            {
              FeaturesWithYLamellae.push_back(gnum);
              size_t distanceFromYmin = j - corners->getComponent(gnum, 1);
              size_t lamellarNumber = distanceFromYmin / 5;
              size_t phaseSelection = lamellarNumber % 2;
              // uint64_t midPoint = corners->getComponent(gnum, 1) + (corners->getComponent(gnum, 4) - corners->getComponent(gnum, 1)) / 2;

              if(phaseSelection == 0)
              {
                m_CellPhases[cellNumber] = 1;
                m_EulerAngles[cellNumber * 3 + 1] = m_EulerAngles[gnum * 3 + 1] + 1.0;
              }
              else
              {
                m_CellPhases[cellNumber] = 2;
                m_EulerAngles[cellNumber * 3 + 1] = m_EulerAngles[gnum * 3 + 1] - 1.0;
                if(m_EulerAngles[cellNumber * 3 + 1] < 0.0)
                  m_EulerAngles[cellNumber * 3 + 1] = 0.0;
              }
            }
            else if(count % 10 == 3) // 3rd, 6th, 9th.... grains with phase 2 will go for Z-Plane-Lamellae
            {
              FeaturesWithZLamellae.push_back(gnum);
              size_t distanceFromZmin = i - corners->getComponent(gnum, 2);
              size_t lamellarNumber = distanceFromZmin / 5;
              size_t phaseSelection = lamellarNumber % 2;
              // uint64_t midPoint = corners->getComponent(gnum, 2) + (corners->getComponent(gnum, 5) - corners->getComponent(gnum, 2)) / 2;

              if(phaseSelection == 0)
              {
                m_CellPhases[cellNumber] = 1;
                m_EulerAngles[cellNumber * 3 + 2] = m_EulerAngles[gnum * 3 + 2] + 1.0;
              }
              else
              {
                m_CellPhases[cellNumber] = 2;
                m_EulerAngles[cellNumber * 3 + 2] = m_EulerAngles[gnum * 3 + 2] - 1.0;
                if(m_EulerAngles[cellNumber * 3 + 2] < 0.0)
                  m_EulerAngles[cellNumber * 3 + 2] = 0.0;
              }
            }
            else if(count % 10 == 4) // 45degree lamellae along XY plane
            {
              FeaturesWithXY45degreeLamellae.push_back(gnum);
              int refXY = (int)k - (int)j;
              if(lamellarDecisionXYplane(refXY) == true)
              {
                m_CellPhases[cellNumber] = 2;
                m_EulerAngles[cellNumber * 3] = m_EulerAngles[gnum * 3] + 1.0;
              }
              else // probably this part is not necessary, as the first cell always considered for the lamellae
              {
                m_CellPhases[cellNumber] = 1;
                m_EulerAngles[cellNumber * 3] = m_EulerAngles[gnum * 3] - 1.0;
                if(m_EulerAngles[cellNumber * 3] < 0.0)
                  m_EulerAngles[cellNumber * 3] = 0.0;
              }
            }
            else if(count % 10 == 5) //-45degree lamellae along XY plane
            {
              FeaturesWithXYminus45degreeLamellae.push_back(gnum);
              int refXY = (int)k + (int)j;
              if(lamellarDecisionXYplane(refXY) == true)
              {
                m_CellPhases[cellNumber] = 2;
                m_EulerAngles[cellNumber * 3] = m_EulerAngles[gnum * 3] + 1.0;
              }
              else // probably this part is not necessary, as the first cell always considered for the lamellae
              {
                m_CellPhases[cellNumber] = 1;
                m_EulerAngles[cellNumber * 3] = m_EulerAngles[gnum * 3] - 1.0;
                if(m_EulerAngles[cellNumber * 3] < 0.0)
                  m_EulerAngles[cellNumber * 3] = 0.0;
              }
            }
            else if(count % 10 == 6) // 45degree lamellae along YZ plane
            {
              FeaturesWithYZ45degreeLamellae.push_back(gnum);
              int refYZ = (int)j - (int)i;
              if(lamellarDecisionXYplane(refYZ) == true)
              {
                m_CellPhases[cellNumber] = 2;
                m_EulerAngles[cellNumber * 3 + 1] = m_EulerAngles[gnum * 3 + 1] + 1.0;
              }
              else // probably this part is not necessary, as the first cell always considered for the lamellae
              {
                m_CellPhases[cellNumber] = 1;
                m_EulerAngles[cellNumber * 3 + 1] = m_EulerAngles[gnum * 3 + 1] - 1.0;
                if(m_EulerAngles[cellNumber * 3 + 1] < 0.0)
                  m_EulerAngles[cellNumber * 3 + 1] = 0.0;
              }
            }
            else if(count % 10 == 7) // minus45degree lamellae along YZ plane
            {
              FeaturesWithYZminus45degreeLamellae.push_back(gnum);
              int refYZ = (int)j + (int)i;
              if(lamellarDecisionXYplane(refYZ) == true)
              {
                m_CellPhases[cellNumber] = 2;
                m_EulerAngles[cellNumber * 3 + 1] = m_EulerAngles[gnum * 3 + 1] + 1.0;
              }
              else // probably this part is not necessary, as the first cell always considered for the lamellae
              {
                m_CellPhases[cellNumber] = 1;
                m_EulerAngles[cellNumber * 3 + 1] = m_EulerAngles[gnum * 3 + 1] - 1.0;
                if(m_EulerAngles[cellNumber * 3 + 1] < 0.0)
                  m_EulerAngles[cellNumber * 3 + 1] = 0.0;
              }
            }
            else if(count % 10 == 8) // 45degree lamellae along ZX plane
            {
              FeaturesWithZX45degreeLamellae.push_back(gnum);
              int refYZ = (int)i - (int)k;
              if(lamellarDecisionXYplane(refYZ) == true)
              {
                m_CellPhases[cellNumber] = 2;
                m_EulerAngles[cellNumber * 3 + 2] = m_EulerAngles[gnum * 3 + 2] + 1.0;
              }
              else // probably this part is not necessary, as the first cell always considered for the lamellae
              {
                m_CellPhases[cellNumber] = 1;
                m_EulerAngles[cellNumber * 3 + 2] = m_EulerAngles[gnum * 3 + 2] - 1.0;
                if(m_EulerAngles[cellNumber * 3 + 2] < 0.0)
                  m_EulerAngles[cellNumber * 3 + 2] = 0.0;
              }
            }
            else // if(count % 10 == 9) // minus45degree lamellae along ZX plane
            {
              FeaturesWithZXminus45degreeLamellae.push_back(gnum);
              int refYZ = (int)i + (int)k;
              if(lamellarDecisionXYplane(refYZ) == true)
              {
                m_CellPhases[cellNumber] = 2;
                m_EulerAngles[cellNumber * 3 + 2] = m_EulerAngles[gnum * 3 + 2] + 1.0;
              }
              else // probably this part is not necessary, as the first cell always considered for the lamellae
              {
                m_CellPhases[cellNumber] = 1;
                m_EulerAngles[cellNumber * 3 + 2] = m_EulerAngles[gnum * 3 + 2] - 1.0;
                if(m_EulerAngles[cellNumber * 3 + 2] < 0.0)
                  m_EulerAngles[cellNumber * 3 + 2] = 0.0;
              }
            }
            /*
            else //XY plane borabor 60degree lamellae
            {
            FeaturesWithXY60degreeLamellae.push_back(gnum);
            RefCellForXY60degreeLamellae.push_back(k);
            RefCellForXY60degreeLamellae.push_back(j);

            m_CellPhases[cellNumber] = 1;
            m_EulerAngles[cellNumber * 3 ] = m_EulerAngles[gnum * 3 ] + 1.0;
          

            }
            */
          }
        }
      }
    }
  }
  //*****************************Generating Lamellae***** end*****************************
}

// -----------------------------------------------------------------------------
// Helper Method - Grabs Index From Matrix Coordinates
// -----------------------------------------------------------------------------
size_t GenerateLamellaeForDuplex::sub2ind(QVector<size_t> tDims, size_t x, size_t y, size_t z) const
{
  return (tDims[1] * tDims[0] * z) + (tDims[0] * y) + x;
}
// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
bool GenerateLamellaeForDuplex::lamellarDecisionXYplane(int x) const
{
  if(x < 0)
    x = -x;
  if(x >= 0 && x <= 2)
    return true;
  else if(((x - 3) / 5) % 2 == 1)
    return true;
  else
    return false;
}
// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void GenerateLamellaeForDuplex::execute()
{
  setErrorCondition(0);
  setWarningCondition(0);
  dataCheck();
  if(getErrorCondition() < 0)
  {
    return;
  }

  generate_lamellaeX();

  notifyStatusMessage(getHumanLabel(), "Complete");
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
AbstractFilter::Pointer GenerateLamellaeForDuplex::newFilterInstance(bool copyFilterParameters) const
{
  GenerateLamellaeForDuplex::Pointer filter = GenerateLamellaeForDuplex::New();
  if(copyFilterParameters)
  {
    copyFilterParameterInstanceVariables(filter.get());
  }
  return filter;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
const QString GenerateLamellaeForDuplex::getCompiledLibraryName() const
{
  return GenericConstants::GenericBaseName;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
const QString GenerateLamellaeForDuplex::getBrandingString() const
{
  return "Generic";
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
const QString GenerateLamellaeForDuplex::getFilterVersion() const
{
  QString version;
  QTextStream vStream(&version);
  vStream << Generic::Version::Major() << "." << Generic::Version::Minor() << "." << Generic::Version::Patch();
  return version;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
const QString GenerateLamellaeForDuplex::getGroupName() const
{
  return SIMPL::FilterGroups::Generic;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
const QString GenerateLamellaeForDuplex::getSubGroupName() const
{
  return SIMPL::FilterSubGroups::MorphologicalFilters;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
const QString GenerateLamellaeForDuplex::getHumanLabel() const
{
  return "Generate Lamellae Along X-Y-Z Planes";
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
const QUuid GenerateLamellaeForDuplex::getUuid()
{
  return QUuid("{d770d728-e291-5de5-8e99-f4cb2e80cc9d}");
}
