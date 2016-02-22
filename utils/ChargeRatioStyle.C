void ChargeRatioStyle()
{
  TStyle* GloStyle;
  GloStyle = gStyle;

  TStyle* ChargeRatioStyle = new TStyle("ChargeRatioStyle", "ChargeRatioStyle");
  gStyle = ChargeRatioStyle;

  //----------------------------------------------------------------------------
  // Canvas
  //----------------------------------------------------------------------------
  ChargeRatioStyle->SetCanvasBorderMode(  0);
  ChargeRatioStyle->SetCanvasBorderSize( 10);
  ChargeRatioStyle->SetCanvasColor     (  0);
  ChargeRatioStyle->SetCanvasDefH      (600);
  ChargeRatioStyle->SetCanvasDefW      (550);
  ChargeRatioStyle->SetCanvasDefX      ( 10);
  ChargeRatioStyle->SetCanvasDefY      ( 10);

  //----------------------------------------------------------------------------
  // Pad
  //----------------------------------------------------------------------------
  ChargeRatioStyle->SetPadBorderMode  (   0);
  ChargeRatioStyle->SetPadBorderSize  (  10);
  ChargeRatioStyle->SetPadColor       (   0);
  ChargeRatioStyle->SetPadBottomMargin(0.20);
  ChargeRatioStyle->SetPadTopMargin   (0.08);
  ChargeRatioStyle->SetPadLeftMargin  (0.18);
  ChargeRatioStyle->SetPadRightMargin (0.05);

  //----------------------------------------------------------------------------
  // Frame
  //----------------------------------------------------------------------------
  ChargeRatioStyle->SetFrameFillStyle ( 0);
  ChargeRatioStyle->SetFrameFillColor ( 0);
  ChargeRatioStyle->SetFrameLineColor ( 1);
  ChargeRatioStyle->SetFrameLineStyle ( 0);
  ChargeRatioStyle->SetFrameLineWidth ( 2);
  ChargeRatioStyle->SetFrameBorderMode( 0);
  ChargeRatioStyle->SetFrameBorderSize(10);

  //----------------------------------------------------------------------------
  // Hist
  //----------------------------------------------------------------------------
  ChargeRatioStyle->SetHistFillColor(0);
  ChargeRatioStyle->SetHistFillStyle(1);
  ChargeRatioStyle->SetHistLineColor(1);
  ChargeRatioStyle->SetHistLineStyle(0);
  ChargeRatioStyle->SetHistLineWidth(1);

  //----------------------------------------------------------------------------
  // Axis
  //----------------------------------------------------------------------------
  ChargeRatioStyle->SetLabelFont  (   42, "xyz");
  ChargeRatioStyle->SetLabelOffset(0.015, "xyz");
  ChargeRatioStyle->SetLabelSize  (0.050, "xyz");
  ChargeRatioStyle->SetNdivisions (  505, "xyz");
  ChargeRatioStyle->SetTitleFont  (   42, "xyz");
  ChargeRatioStyle->SetTitleSize  (0.050, "xyz");

  ChargeRatioStyle->SetTitleOffset(  1.4,   "x");
  ChargeRatioStyle->SetTitleOffset(  1.2,   "y");
  ChargeRatioStyle->SetPadTickX   (           1);  // Tick marks on the opposite side of the frame
  ChargeRatioStyle->SetPadTickY   (           1);  // Tick marks on the opposite side of the frame


  //----------------------------------------------------------------------------
  // Title
  //----------------------------------------------------------------------------
  ChargeRatioStyle->SetTitleBorderSize( 0);
  ChargeRatioStyle->SetTitleFillColor (10);
  ChargeRatioStyle->SetTitleFont      (42, "");

  //----------------------------------------------------------------------------
  // Stat
  //----------------------------------------------------------------------------
  ChargeRatioStyle->SetOptStat       (1110);
  ChargeRatioStyle->SetStatBorderSize(   0);
  ChargeRatioStyle->SetStatColor     (  10);
  ChargeRatioStyle->SetStatFont      (  42);
  ChargeRatioStyle->SetStatX         (0.94);
  ChargeRatioStyle->SetStatY         (0.91);

  return();
}
