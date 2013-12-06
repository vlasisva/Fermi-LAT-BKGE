//Author: Vlasios Vasileiou <vlasisva@gmail.com>
void prettyroot(){
 gStyle->SetPalette(1);
 gStyle->SetCanvasColor(kWhite);     // background is no longer mouse-dropping white. 
 gStyle->SetOptStat(0);

 gStyle->SetPadTopMargin(0.07);
 gStyle->SetPadLeftMargin(0.13);
 gStyle->SetPadRightMargin(0.11);
 gStyle->SetPadBottomMargin(0.1);
 gStyle->SetPadTickX(1);  //make ticks be on all 4 sides.
 gStyle->SetPadTickY(1);
 gStyle->SetOptTitle(1);


  gStyle->SetTitleFont(42,"xyz");
  gStyle->SetTitleSize(0.08,"xy");
  gStyle->SetTitleFillColor(0);
//  gStyle->SetTitleOffset(1.4,"x");
//  gStyle->SetTitleOffset(1,"y");

 // statistics box
  gStyle->SetStatFont(42);
  gStyle->SetStatX(.91);
  gStyle->SetStatY(.90);
  gStyle->SetStatW(.15);
  gStyle->SetStatH(.15);
  gStyle->SetGridColor(16);

  gStyle->SetLegendBorderSize(1);
  gStyle->SetFrameBorderMode(0);
  gStyle->SetLineScalePS(1);
  gStyle->SetCanvasBorderMode(0);     // turn off canvas borders. 
  gStyle->SetPadBorderMode(0);
  gStyle->SetPaintTextFormat("5.2f");  // What precision to put numbers if plotted with "TEXT". 

  // For publishing:
  gStyle->SetLineWidth(1.5);
//  gStyle->SetTextSize(0.9);
  gStyle->SetLabelSize(0.05,"xy");
  gStyle->SetLabelFont(42,"xyz");
  gStyle->SetLabelOffset(0.01,"xy");
  gStyle->SetPadTopMargin(0.1);
  gStyle->SetPadRightMargin(0.1);
  gStyle->SetPadBottomMargin(0.16);
  gStyle->SetPadLeftMargin(0.12);


}
