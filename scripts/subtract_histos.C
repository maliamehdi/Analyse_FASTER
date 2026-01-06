// subtract_histos.C
// ROOT macro to normalize two histograms by their integrals and subtract them.
// Usage (ROOT):
//   root -l -b -q 'scripts/subtract_histos.C("fileA.root","histA","fileB.root","histB","out.root","h_sub", true, 2, 1.0)'
// - fileA.root, fileB.root: input ROOT files (can include directories in hist names like "dir1/hist")
// - histA, histB: histogram names/paths inside the files
// - out.root: output filename (created/recreated)
// - h_sub: name for the output subtracted histogram (optional, defaults to "h_sub")
// - saveInputs: also store the normalized inputs into the output file (default true)
// - scaleWhich: apply scaleFactor to 1 = histA, 2 = histB (default 2)
// - scaleFactor: multiplicative factor applied AFTER normalization to the chosen histogram (default 1.0)
// Notes:
// - Histograms must be 1D. The macro will automatically:
//   1) Align ranges by restricting to the maximal common range supported by both histograms
//      (i.e. use common bin edges; effectively cropping the larger to the smaller when possible)
//   2) Rebin minimally so both have identical bin edges using the intersection of their edge sets
//      (works for variable-bin histograms and uniform ones; preserves errors)
//   3) Normalize each by its integral (width-weighted), then apply `scaleFactor` to only one histogram,
//      and finally subtract: result = normA - (scaleFactor * normB) if scaleWhich=2 (default)
// - If you truly want to scale BEFORE normalization, note that it cancels out during normalization; instead
//   use `scaleFactor` as provided (applied after normalization) to get a meaningful effect.

#include <TFile.h>
#include <TH1.h>
#include <TKey.h>
#include <TClass.h>
#include <TROOT.h>
#include <TString.h>
#include <TError.h>

#include <iostream>
#include <cmath>
#include <memory>

namespace {
  double safeIntegral(const TH1* h) {
    if (!h) return 0.0;
    // Prefer width-weighted integral for variable bin width robustness
    double integ = h->Integral("width");
#include <vector>
#include <algorithm>
    if (integ <= 0.0 || !std::isfinite(integ)) {
      // Fallback to simple sum of bin contents excluding under/overflow
  constexpr double kEps = 1e-9;

  inline bool almostEqual(double a, double b, double eps = kEps) {
    double scale = std::max(1.0, std::max(std::fabs(a), std::fabs(b)));
    return std::fabs(a - b) <= eps * scale;
  }

      integ = h->Integral(1, h->GetNbinsX());
    }
    return integ;
  }

  bool sameBinning(const TH1* a, const TH1* b) {
    if (!a || !b) return false;
    if (a->GetDimension() != 1 || b->GetDimension() != 1) return false;
    const int na = a->GetNbinsX();
    const int nb = b->GetNbinsX();
    if (na != nb) return false;
  bool sameBinning(const TH1* a, const TH1* b) {
    if (!a || !b) return false;
    if (a->GetDimension() != 1 || b->GetDimension() != 1) return false;
    const int na = a->GetNbinsX();
    const int nb = b->GetNbinsX();
    if (na != nb) return false;
    const TAxis* xa = a->GetXaxis();
    const TAxis* xb = b->GetXaxis();
    // Compare all bin edges (including upper edge of last bin)
    for (int i = 1; i <= na; ++i) {
      double la = xa->GetBinLowEdge(i);
      double lb = xb->GetBinLowEdge(i);
      if (!almostEqual(la, lb)) return false;
    }
    double ua = xa->GetBinUpEdge(na);
    double ub = xb->GetBinUpEdge(nb);
    if (!almostEqual(ua, ub)) return false;
    return true;
  }

  // Extract full set of x bin edges for a histogram (size nbins+1)
  std::vector<double> getEdges(const TH1* h) {
    std::vector<double> edges;
    if (!h) return edges;
    const TAxis* ax = h->GetXaxis();
    const int nb = ax->GetNbins();
    const TArrayD* xb = ax->GetXbins();
    if (xb && xb->GetSize() > 0) {
      edges.reserve(xb->GetSize());
      for (int i = 0; i < xb->GetSize(); ++i) edges.push_back((*xb)[i]);
    } else {
      // Uniform bins
      edges.reserve(nb + 1);
      const double xmin = ax->GetXmin();
      const double xmax = ax->GetXmax();
      const double step = (xmax - xmin) / nb;
      for (int i = 0; i <= nb; ++i) edges.push_back(xmin + i * step);
    }
    return edges;
  }

  // Compute intersection of two sorted edge arrays with tolerance
  std::vector<double> intersectEdges(const std::vector<double>& a, const std::vector<double>& b) {
    std::vector<double> out;
    size_t i = 0, j = 0;
    out.reserve(std::min(a.size(), b.size()));
    while (i < a.size() && j < b.size()) {
      double va = a[i], vb = b[j];
      if (almostEqual(va, vb)) { out.push_back(va); ++i; ++j; }
      else if (va < vb) { ++i; }
      else { ++j; }
    }
    // Ensure strictly increasing and unique (guard against duplicates due to tolerance)
    std::vector<double> unique;
    unique.reserve(out.size());
    for (double v : out) {
      if (unique.empty() || !almostEqual(unique.back(), v)) unique.push_back(v);
    }
    return unique;
  }

  // Rebin a histogram to the provided edges using TH1::Rebin(xbins) API
  std::unique_ptr<TH1> rebinToEdges(TH1* h, const std::vector<double>& edges, const char* newname) {
    if (!h) return nullptr;
    if (edges.size() < 2) return nullptr;
    // TH1::Rebin with xbins expects number of new bins and pointer to edges array
    // Provide a new name so the original is left unchanged by ROOT
    std::unique_ptr<TH1> result( h->Rebin(static_cast<Int_t>(edges.size()-1), newname, edges.data()) );
    if (result) result->SetDirectory(nullptr);
    return result;
  }
    return h;
  }
}

int subtract_histos(const char* fileA,
                    const char* histA,
                    const char* fileB,
                    const char* outHistName = "h_sub",
                    bool saveInputs = true,
                    int scaleWhich = 2,
                    double scaleFactor = 1.0,
                    const char* outHistName = "h_sub",
                    bool saveInputs = true) {
  if (!fileA || !histA || !fileB || !histB || !outFile) {
    ::Error("subtract_histos", "Invalid arguments (null string)");
    return 1;
  }

  std::unique_ptr<TFile> fa(TFile::Open(fileA, "READ"));
  if (!fa || fa->IsZombie()) {
    ::Error("subtract_histos", "Cannot open input file A: %s", fileA);
    return 1;
  }

  std::unique_ptr<TFile> fb(TFile::Open(fileB, "READ"));
  if (!fb || fb->IsZombie()) {
    ::Error("subtract_histos", "Cannot open input file B: %s", fileB);
    return 1;
  }

  TH1* hA = fetchHist(fa.get(), histA);
  if (!hA) {
    ::Error("subtract_histos", "Histogram not found in A: %s", histA);
    return 1;
  }
  TH1* hB = fetchHist(fb.get(), histB);
  if (!hB) {
    ::Error("subtract_histos", "Histogram not found in B: %s", histB);
    return 1;
  }

  if (hA->GetDimension() != 1 || hB->GetDimension() != 1) {
    ::Error("subtract_histos", "Only 1D histograms are supported");
    return 1;
  }

  // Clone inputs to detached working copies
  std::unique_ptr<TH1> hA_work(static_cast<TH1*>(hA->Clone("hA_work"))); hA_work->SetDirectory(nullptr);
  std::unique_ptr<TH1> hB_work(static_cast<TH1*>(hB->Clone("hB_work"))); hB_work->SetDirectory(nullptr);

  // Build full edge arrays
  auto edgesA = getEdges(hA_work.get());
  auto edgesB = getEdges(hB_work.get());
  // Compute common edges (this simultaneously aligns ranges and prepares minimal common binning)
  auto commonEdges = intersectEdges(edgesA, edgesB);
  if (commonEdges.size() < 2) {
    ::Error("subtract_histos", "No common bin edges found between histograms; cannot align range/binning.");
    return 2;
  }

  // Rebin both to the common edges (this also crops to the common range)
  auto hA_aligned = rebinToEdges(hA_work.get(), commonEdges, "hA_aligned");
  auto hB_aligned = rebinToEdges(hB_work.get(), commonEdges, "hB_aligned");
  if (!hA_aligned || !hB_aligned) {
    ::Error("subtract_histos", "Failed to rebin to common edges.");
    return 2;
  }

  // Normalize by integrals (width-weighted)
  const double iA = safeIntegral(hA_aligned.get());
  const double iB = safeIntegral(hB_aligned.get());
  if (iA <= 0.0 || !std::isfinite(iA)) {
    ::Error("subtract_histos", "Non-positive or invalid integral for histA (%s).", histA);
    return 2;
  }
  if (iB <= 0.0 || !std::isfinite(iB)) {
    ::Error("subtract_histos", "Non-positive or invalid integral for histB (%s).", histB);
    return 2;
  }

  std::unique_ptr<TH1> hA_norm(static_cast<TH1*>(hA_aligned->Clone("hA_norm_temp"))); hA_norm->SetDirectory(nullptr);
  std::unique_ptr<TH1> hB_norm(static_cast<TH1*>(hB_aligned->Clone("hB_norm_temp"))); hB_norm->SetDirectory(nullptr);
  hA_norm->Scale(1.0 / iA);
  hB_norm->Scale(1.0 / iB);

  // Apply additional scaling factor to only one of the histograms AFTER normalization
  if (scaleWhich == 1) hA_norm->Scale(scaleFactor);
  else if (scaleWhich == 2) hB_norm->Scale(scaleFactor);

  // Create the subtracted histogram
  std::unique_ptr<TH1> hSub(static_cast<TH1*>(hA_norm->Clone(outHistName && outHistName[0] ? outHistName : "h_sub")));
  hSub->SetDirectory(nullptr);
  hSub->Add(hB_norm.get(), -1.0);

  // Title with provenance
  TString title;
  title.Form("(%s:%s)/Int - (%s:%s)/Int", fileA, histA, fileB, histB);
  hSub->SetTitle(title);

  // Write to output
  std::unique_ptr<TFile> fout(TFile::Open(outFile, "RECREATE"));
  if (!fout || fout->IsZombie()) {
    ::Error("subtract_histos", "Cannot create output file: %s", outFile);
    return 1;
  }

  hSub->Write();
  if (saveInputs) {
    hA_norm->SetName("hA_norm");
    hB_norm->SetName("hB_norm");
    hA_norm->Write();
    hB_norm->Write();
  }
  fout->Write();
  fout->Close();

  ::Info("subtract_histos", "Wrote subtracted histogram '%s' to %s", hSub->GetName(), outFile);
  return 0;
}
