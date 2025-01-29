#ifndef TEMP_HEADER
#define TEMP_HEADER

namespace oomph
{
  class Temp
  {
  private:
    /// Storage for the doc info
    DocInfo Doc_info;

  public:
    /// Constructor
    Temp() : Doc_info("RESLT") {}

    /// access the doc info
    DocInfo& doc_info()
    {
      return Doc_info;
    }
  };
}; // namespace oomph

#endif
