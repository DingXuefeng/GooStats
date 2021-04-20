/*****************************************************************************/
// Author: Xuefeng Ding <xuefeng.ding.physics@gmail.com>
// Insitute: Gran Sasso Science Institute, L'Aquila, 67100, Italy
// Date: 2018 April 7th
// Version: v1.0
// Description: GooStats, a statistical analysis toolkit that runs on GPU.
//
// All rights reserved. 2018 copyrighted.
/*****************************************************************************/
#ifndef PlotManager_H
#define PlotManager_H
#include <map>
#include <string>
#include <vector>
#include <set>
#include "Rtypes.h"
#include "StatModule.h"
class TObject;
class DatasetManager;
class GooPdf;
class TH1;
class TF1;
class SumPdf;
class TCanvas;
class TFile;
class GSFitManager;
class PlotManager : public StatModule {
  public:
    PlotManager(const std::string &_name="PlotManager") : StatModule(_name) { }
    virtual bool init() override;
    virtual bool finish() override;

  protected:
    struct Config {
      int color; ///< color
      int style; ///< line style , 1 is solid line
      int width; ///< line width
    };
    class TF1Helper {
      public:
        TF1Helper(GooPdf *pdf_,double norm_ = 1,int index = 0);
        TF1 *getTF1() { return f; }
      private:
        double eval(double *xx,double *par);
        TH1 *data = nullptr; // root will delete it
      private:
        TF1 *f = nullptr; // root will delete it
    };
  public:
    virtual void draw(int event,const std::vector<DatasetManager*> &datasets);
    void drawLikelihoodpValue(int event,double LL,const std::vector<double> &LLs);
  public:
    virtual TCanvas *drawSingleGroup(const std::string &name,const std::vector<DatasetManager*> &datasets) ;
    //! user can set the color and line style of each components by specifying the names
    //! plot() will look up the name of components of SumPdf
    //! by default different colors will be used for each species and all solid lines
    virtual void draw(GSFitManager *gsFitManager,SumPdf *pdf,std::map<std::string,Config> config= std::map<std::string,Config>(),int index = 0);
  protected:
    EColor getColor(const std::string &n) const {
      return colorlibrary.find(n)!=colorlibrary.end()?colorlibrary.at(n):kBlack;
    }
    bool createPdf() const { return _createPdf; };
    const std::string &outName() const;
    std::set<TObject*> toBeSaved;
    virtual void set_gStyle();
  private:
    bool _createPdf = true;
    //! configuration for plooting all species
    std::map<std::string, EColor> colorlibrary {
      {"kRed", kRed},
        {"kOrange", kOrange},
        {"kYellow", kYellow},
        {"kGreen", kGreen},
        {"kBlue", kBlue},
        {"kViolet", kViolet},
        {"kAzure", kAzure},
        {"kWhite", kWhite},
        {"kBlack", kBlack},
        {"kGray", kGray},
        {"kMagenta", kMagenta},
        {"kCyan", kCyan},
        {"kSpring", kSpring},
        {"kTeal", kTeal},
        {"kPink", kPink},
        {"kBlack", kBlack},
    };
};
#endif
