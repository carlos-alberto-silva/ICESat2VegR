#include <vector>
#include <Rcpp.h>
#include <gdal_priv.h>
#include <proj.h>
#include "cpl_string.h"
#include <string>

using namespace Rcpp;

// Function to get the file extension from a file name
std::string GetFileExtension(std::string fileName)
{
  size_t dotPos = fileName.find_last_of('.');
  if (dotPos == std::string::npos)
  {
    return "";
  }
  return fileName.substr(dotPos + 1);
}

GDALDriver *GetDriverByFile(std::string filename)
{
  // Get the file extension from the command line argument
  std::string fileExtension = GetFileExtension(filename);

  if (fileExtension.empty())
  {
    Rcpp::stop("Invalid file name or no extension found: %s", filename);
  }

  GDALDriver *pDriver = nullptr;
  GDALDriverManager *driverManager = GetGDALDriverManager();
  int nDrivers = driverManager->GetDriverCount();
  for (int i = 0; i < nDrivers; i++)
  {
    GDALDriver *tempDriver = driverManager->GetDriver(i);
    const char *driverExtensions = tempDriver->GetMetadataItem(GDAL_DMD_EXTENSIONS);

    if (driverExtensions != nullptr)
    {
      // Split the extensions string and compare with the given file extension
      char **tokens = CSLTokenizeString(driverExtensions);
      for (int j = 0; tokens[j] != nullptr; j++)
      {
        if (EQUAL(tokens[j], fileExtension.c_str()))
        {
          pDriver = tempDriver;
          break;
        }
      }
      CSLDestroy(tokens);
    }

    if (pDriver != nullptr)
    {
      break;
    }
  }

  const char *rasterCapability = pDriver->GetMetadataItem(GDAL_DCAP_RASTER);

  if (!(rasterCapability && EQUAL(rasterCapability, "YES")))
  {
    Rcpp::stop("Driver %s does not support raster data", pDriver->GetDescription());
  }

  if (pDriver == nullptr)
  {
    Rcpp::stop("No driver found for extension '%s'", fileExtension);
  }

  return pDriver;
}

class GDALRasterBandR
{
private:
  GDALRasterBand *band;

public:
  GDALRasterBandR(GDALRasterBand *the_band)
  {
    band = the_band;
  }

  void CalculateStatistics()
  {
    band->ComputeStatistics(false, NULL, NULL, NULL, NULL, NULL, NULL);
  }

  IntegerVector GetBlockXSize()
  {
    IntegerVector result(1);
    band->GetBlockSize(result.begin(), NULL);
    return result;
  }

  IntegerVector GetBlockYSize()
  {
    IntegerVector result(1);
    band->GetBlockSize(NULL, result.begin());
    return result;
  }

  int GetRasterDataType()
  {
    return (int)band->GetRasterDataType();
  }

  int GetXSize()
  {
    return band->GetXSize();
  }

  int GetYSize()
  {
    return band->GetYSize();
  }

  double GetNoDataValue()
  {
    return band->GetNoDataValue();
  }

  template <typename T, typename S>
  S ReadBlock(int iXBlock, int iYBlock)
  {
    CPLErr res = CE_None;
    S output;

    int nXBlockSize, nYBlockSize;
    band->GetBlockSize(&nXBlockSize, &nYBlockSize);
    if (std::is_same<T, GInt32>::value || std::is_same<T, double>::value)
    {
      S vec(nXBlockSize * nYBlockSize);
      res = band->ReadBlock(iXBlock, iYBlock, vec.begin());
      output = vec;
    }
    else
    {
      std::vector<T> buffer(nXBlockSize * nYBlockSize);
      res = band->ReadBlock(iXBlock, iYBlock, buffer.data());
      output = Rcpp::wrap(buffer.begin(), buffer.end());
    }

    if (res == CE_Failure)
      Rcpp::stop(CPLGetLastErrorMsg());

    return output;
  }

  template <typename T>
  void WriteBlock(int iXBlock, int iYBlock, T buffer)
  {
    GDALDataType dtype = band->GetRasterDataType();
    CPLErr res = CE_None;

    switch (dtype)
    {
    case GDALDataType::GDT_UInt16:
    {
      std::vector<GUInt16> vec(buffer.begin(), buffer.end());
      res = band->WriteBlock(iXBlock, iYBlock, vec.data());
      break;
    }

    case GDALDataType::GDT_Int16:
    {
      std::vector<GInt16> vec(buffer.begin(), buffer.end());
      res = band->WriteBlock(iXBlock, iYBlock, vec.data());
      break;
    }

    case GDALDataType::GDT_UInt32:
    {
      std::vector<GUInt32> vec(buffer.begin(), buffer.end());
      res = band->WriteBlock(iXBlock, iYBlock, vec.data());
      break;
    }

    case GDALDataType::GDT_Float32:
    {
      std::vector<float> vec(buffer.begin(), buffer.end());
      res = band->WriteBlock(iXBlock, iYBlock, vec.data());
      break;
    }

    default:
    {
      res = band->WriteBlock(iXBlock, iYBlock, buffer.begin());
      break;
    }
    }

    if (res == CE_Failure)
      Rcpp::stop(CPLGetLastErrorMsg());
  }
};

class GDALDatasetR
{
private:
  GDALDataset *ds = NULL;
  bool closed = false;

public:
  GDALDatasetR(GDALDataset *_ds)
  {
    ds = _ds;
  }

  virtual ~GDALDatasetR()
  {
    ds = NULL;
  }

  GDALRasterBandR *GetRasterBand(int nband)
  {
    GDALRasterBandR *band = new GDALRasterBandR(ds->GetRasterBand(nband));
    return band;
  }

  int GetRasterXSize()
  {
    return ds->GetRasterXSize();
  }

  int GetRasterYSize()
  {
    return ds->GetRasterYSize();
  }

  void Close()
  {
    if (!closed)
    {
      GDALClose((GDALDatasetH)ds);
      closed = true;
    }
  }
};

GDALDatasetR *create_dataset(
    std::string output,
    int nbands,
    int datatype,
    std::string projection,
    double lat_min,
    double lat_max,
    double lon_min,
    double lon_max,
    std::vector<double> res,
    double nodata,
    CharacterVector co)
{
  CPLErr err = CE_None;
  int width = (int)ceil((lon_max - lon_min) / res[0]);
  int height = (int)ceil((lat_min - lat_max) / res[1]);

  GDALDriver *driver = GetDriverByFile(output);
  if (driver == nullptr)
    Rcpp::stop(CPLGetLastErrorMsg());

  std::vector<char *> charVec{};

  for (auto &option : co)
  {
    charVec.push_back(option);
  }
  charVec.push_back(nullptr);

  GDALDataset *ds = driver->Create(output.c_str(), width, height, nbands, (GDALDataType)datatype, charVec.data());

  if (ds == NULL)
    Rcpp::stop(CPLGetLastErrorMsg());

  double transform[6] = {lon_min, res[0], 0, lat_max, 0, res[1]};
  ds->SetGeoTransform(transform);
  for (int i = 1; i <= nbands; i++)
  {
    GDALRasterBand *band = ds->GetRasterBand(i);
    err = band->SetNoDataValue(nodata);
    if (err == CE_Failure)
      Rcpp::stop(CPLGetLastErrorMsg());
  }
  err = ds->SetProjection(projection.c_str());

  if (err == CE_Failure)
    Rcpp::stop(CPLGetLastErrorMsg());

  GDALDatasetR *outDs = new GDALDatasetR(ds);

  return outDs;
}

void GDALDatasetFinalizer(GDALDatasetR *ds)
{
  ds->Close();
}


void InitializeGDAL(std::vector<std::string> paths)
{
  std::vector<const char*> cstrings;
  cstrings.reserve(paths.size());

  for (auto &str : paths)
  {
    cstrings.push_back(str.c_str());
  }

  proj_context_set_search_paths(NULL, cstrings.size(), cstrings.data());

  GDALAllRegister();
  CPLSetErrorHandler(CPLQuietErrorHandler);
}

IntegerVector GetProjVersion()
{
  PJ_INFO info = proj_info();
  int major, minor, patch;
  if (sscanf(info.version, "%d.%d.%d", &major, &minor, &patch) != 3)
  {
    Rcpp::stop("Failed to parse PROJ version");
  }
  return IntegerVector::create(major, minor, patch);
}

GDALDatasetR *RGDALOpen(const char *filename, bool readonly)
{
  GDALAccess accessMode = readonly ? GDALAccess::GA_ReadOnly : GDALAccess::GA_Update;
  GDALDataset *ds = (GDALDataset *)GDALOpen(filename, accessMode);
  if (ds == NULL)
  {
    Rcpp::stop(CPLGetLastErrorMsg());
  }
  GDALDatasetR *outDs = new GDALDatasetR(ds);

  return outDs;
}

RCPP_MODULE(gdal_module)
{
  class_<GDALDatasetR>("CPP_GDALDataset")
      .method("GetRasterBand", &GDALDatasetR::GetRasterBand)
      .method("GetRasterXSize", &GDALDatasetR::GetRasterXSize)
      .method("GetRasterYSize", &GDALDatasetR::GetRasterYSize)
      .method("Close", &GDALDatasetR::Close)
      .finalizer(&GDALDatasetFinalizer);

  class_<GDALRasterBandR>("CPP_GDALRasterBand")
      .method("ReadBlock1", &GDALRasterBandR::ReadBlock<GByte, RawVector>)
      .method("ReadBlock2", &GDALRasterBandR::ReadBlock<GUInt16, IntegerVector>)
      .method("ReadBlock3", &GDALRasterBandR::ReadBlock<GInt16, IntegerVector>)
      .method("ReadBlock4", &GDALRasterBandR::ReadBlock<GUInt32, IntegerVector>)
      .method("ReadBlock5", &GDALRasterBandR::ReadBlock<GInt32, IntegerVector>)
      .method("ReadBlock6", &GDALRasterBandR::ReadBlock<float, NumericVector>)
      .method("ReadBlock7", &GDALRasterBandR::ReadBlock<double, NumericVector>)
      .method("WriteBlock1", &GDALRasterBandR::WriteBlock<RawVector>)
      .method("WriteBlock2", &GDALRasterBandR::WriteBlock<IntegerVector>)
      .method("WriteBlock3", &GDALRasterBandR::WriteBlock<IntegerVector>)
      .method("WriteBlock4", &GDALRasterBandR::WriteBlock<IntegerVector>)
      .method("WriteBlock5", &GDALRasterBandR::WriteBlock<IntegerVector>)
      .method("WriteBlock6", &GDALRasterBandR::WriteBlock<NumericVector>)
      .method("WriteBlock7", &GDALRasterBandR::WriteBlock<NumericVector>)
      .method("GetBlockXSize", &GDALRasterBandR::GetBlockXSize)
      .method("GetBlockYSize", &GDALRasterBandR::GetBlockYSize)
      .method("GetNoDataValue", &GDALRasterBandR::GetNoDataValue)
      .method("GetXSize", &GDALRasterBandR::GetXSize)
      .method("GetYSize", &GDALRasterBandR::GetYSize)
      .method("GetRasterDataType", &GDALRasterBandR::GetRasterDataType)
      .method("CalculateStatistics", &GDALRasterBandR::CalculateStatistics);

  function("create_dataset", &create_dataset);
  function("RGDALOpen", &RGDALOpen);
  function("InitializeGDAL", &InitializeGDAL);
  function("GetProjVersion", &GetProjVersion);
}
