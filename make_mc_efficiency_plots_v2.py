import os
import json
import math
import re
import yaml
import numpy as np
import pickle as pkl
from tqdm import tqdm
from pprint import pprint
from pathlib import Path
from uncertainties import ufloat
import argparse
import copy

import ROOT

from plotting_scripts.EfficiencyPlot import EfficiencyPlot


def check_parameter_limits(savename, fit_result):
    if not fit_result:
        return 

    params_at_limit = []
    final_params = fit_result.floatParsFinal()

    for i in range(final_params.getSize()):
        param = final_params.at(i)
        name = param.GetName()
        value = param.getVal()
        min_val = param.getMin()
        max_val = param.getMax()

        # Use a small tolerance for floating point comparison
        tolerance = 1e-6

        if abs(value - min_val) < tolerance:
            params_at_limit.append(f"  - {name} (at MIN limit: {min_val})")
        
        if abs(value - max_val) < tolerance:
            params_at_limit.append(f"  - {name} (at MAX limit: {max_val})")

    if params_at_limit:
        print(f"--- Boundary Check for: {savename} ---")
        print("WARNING: The following parameters have final values at their limits:")
        for param_info in params_at_limit:
            print(param_info)
        print(50*'-')


def set_verbosity(verb):
    ROOT.gROOT.SetBatch(True)
    ROOT.gErrorIgnoreLevel = ROOT.kInfo if verb else ROOT.kError
    ROOT.RooMsgService.instance().setGlobalKillBelow(ROOT.RooFit.INFO if verb else ROOT.RooFit.ERROR)
    printlevel = ROOT.RooFit.PrintLevel(1 if verb else -1)

    return printlevel


def get_event_count(h, xrange=None):
    nEvents = 0.0
    err = 0.0
    
    if xrange is not None:
        bin_low = h.GetXaxis().FindBin(xrange[0])
        bin_high = h.GetXaxis().FindBin(xrange[1])
        
        error_ref = ROOT.Double(0.0)
        
        nEvents = h.IntegralAndError(bin_low, bin_high, error_ref)
        err = float(error_ref)
    
    else:
        nEvents = h.GetSumOfWeights() 
        err = np.sqrt(nEvents)
        
    return nEvents, err


def estimateBins(h,nbins=5):
    evts = [h.GetBinContent(i) for i in range(h.GetXaxis().GetNbins())]
    tot = np.sum(evts)

    frac = 1
    n_evts = 0
    bins = f'{h.GetBinLowEdge(1)}'
    for i in range(1, h.GetXaxis().GetNbins()):
        n_evts += h.GetBinContent(i)
        if n_evts > frac * tot / nbins:
            bins += ', '+str(round(h.GetBinLowEdge(i),2))
            frac += 1
    bins += ', '+str(round(h.GetBinLowEdge(h.GetXaxis().GetNbins()+1),2))
    return bins

def assign_hist_format(h_name):
    if 'pt1pt2binned' in h_name:
        h_format = {
            'bins_x' : np.array([4, 6, 9, 10, 12, 20], dtype=np.double),
            'bins_y' : np.array([4, 6, 9, 10, 12, 20], dtype=np.double),            
            'xlabel' : 'Leading Electron p_{T} [GeV]',
            'ylabel' : 'Subleading Electron p_{T} [GeV]',
            'is_2d' : True,
        }
    elif 'ptbinned' in h_name:
        h_format = {
            'bins' : np.array([4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 20], dtype=np.double),
            # 'bins' : np.array([5, 6, 7, 8, 9, 10, 11, 12, 13, 20],  dtype=np.double),
            # 'bins' : np.array([4, 6.5, 9, 12, 20], dtype=np.double),
            # 'bins' : np.array([4, 6, 9, 10, 12, 20], dtype=np.double),
            #[5, 11, 999
            # 'bins' : np.array([5, 8, 20], dtype=np.double),
            # 'bins' : np.array([5, 8, 12, 20], dtype=np.double),
            'xlabel' : 'Sublead Electron p_{T} [GeV]',
            'is_2d' : False,
        }
    elif 'drbinned' in h_name:
        h_format = {
            'bins' : np.array([0, 0.12, 0.2, 0.28, 0.44, 1.], dtype=np.double),
            'xlabel' : r'\Delta R(e_{1},e_{2})',
            'is_2d' : False,
        }
    elif 'etabinned' in h_name:
        h_format = {
            'bins' : np.array([-1.22, -0.7, -.2, 0.2, .7, 1.22], dtype=np.double),
            'xlabel' : r'Sublead \ Electron \ \eta',
            'is_2d' : False,
        }
    else:
        raise ValueError(f'No formatting options for hist name {h_name}')

    return h_format


def do_fit(h, signal_only=False, savename=None, get_params=False, signal_params=None, bkg_shape_deg=4, printlevel=ROOT.RooFit.PrintLevel(-1), chi2_threshold=15.0, max_tries=20):

    data_yield = h.GetEntries()
    data_min = h.GetXaxis().GetXmin()
    data_max = h.GetXaxis().GetXmax()

    mass = ROOT.RooRealVar('mass', 'mass', 2.5, 3.4)
    data = ROOT.RooDataHist('data', 'data', mass, h)

    sig_coeff =  ROOT.RooRealVar('sig_coeff', 'sig_coeff', data_yield*.03, 0, data_yield)
    cb_mean =   ROOT.RooRealVar('cb_mean', 'cb_mean', 3.0969)
    cb_sigma =  ROOT.RooRealVar('cb_sigma', 'cb_sigma', 0.046, 0.02, 0.1)
    cb_alpha = ROOT.RooRealVar('cb_alpha', 'cb_alpha', 2.0, 0.3, 10.)
    cb_n =     ROOT.RooRealVar('cb_n', 'cb_n', 8.0, 1.01, 40.)
    sig_pdf = ROOT.RooCrystalBall('sig_pdf', 'Signal Fit', mass, cb_mean, cb_sigma, cb_alpha, cb_n)

    bkg_coeff =  ROOT.RooRealVar('bkg_coeff', 'bkg_coeff', data_yield*.97, 0, data_yield)
    bern_pars = [ROOT.RooRealVar(f'bern_p{i}', f'bern_p{i}', 1, 0, 10) for i in range(bkg_shape_deg)]
    bkg_pdf = ROOT.RooBernstein('bkg_pdf', 'Background Fit', mass, ROOT.RooArgList(*bern_pars))

    # Tail constraint
    n_central_val = ROOT.RooConstVar("n_central_val", "n_central_val", 10.0) # Target a large n
    n_sigma_val = ROOT.RooConstVar("n_sigma_val", "n_sigma_val", 3.0)   # With a reasonably strong pull
    n_constraint_pdf = ROOT.RooGaussian(
        "n_constraint_pdf", "Constraint on cb_n",
        cb_n, n_central_val, n_sigma_val
    )

    if signal_only:
        fit_model = ROOT.RooAddPdf('fit_model', 'Signal Fit', ROOT.RooArgList(sig_pdf), ROOT.RooArgList(sig_coeff))
    else:
        fit_model = ROOT.RooAddPdf('fit_model', 'Signal + Background Fit', ROOT.RooArgList(sig_pdf, bkg_pdf), ROOT.RooArgList(sig_coeff, bkg_coeff))


    fit_result = None
    chi2 = float('inf')
    attempt = 0
    sig_unc_ratio = 1
    good_fit = False

    while not good_fit and attempt < max_tries:
        if fit_result:
            fit_result.Delete()

        fit_result = fit_model.fitTo(
            data,
            ROOT.RooFit.Save(True),
            ROOT.RooFit.Minos(True),
            ROOT.RooFit.ExternalConstraints(ROOT.RooArgSet(n_constraint_pdf)),
            printlevel
        )

        fit_status = fit_result.status()
        cov_qual = fit_result.covQual()
        if fit_status >= 0 and cov_qual >= 2:
            sig_unc_ratio = sig_coeff.getError() / sig_coeff.getVal() if sig_coeff.getVal() > 0 else float('inf')
            
            tmp_frame = mass.frame(ROOT.RooFit.Title('tmp'))
            data.plotOn(tmp_frame)
            if signal_only:
                fit_model.plotOn(tmp_frame, ROOT.RooFit.Normalization(sig_coeff.getVal(),ROOT.RooAbsReal.NumEvent))
            else:
                fit_model.plotOn(tmp_frame, ROOT.RooFit.Normalization(sig_coeff.getVal()+bkg_coeff.getVal(),ROOT.RooAbsReal.NumEvent))
            
            chi2 = tmp_frame.chiSquare('fit_model_Norm[mass]', 'h_data', len(fit_result.floatParsFinal()))
            tmp_frame.Delete()

            if not math.isnan(chi2) and chi2 < chi2_threshold and sig_unc_ratio < 1:
                good_fit = True
            
        if fit_result and not good_fit:
            params_to_randomize = [cb_sigma, cb_alpha, cb_n] + bern_pars 
            for param in params_to_randomize:
                current_val = param.getVal()
                min_val = param.getMin()
                max_val = param.getMax()
                
                perturbation_range = (max_val - min_val) * 0.1
                random_shift = np.random.uniform(-perturbation_range, perturbation_range)
                new_val = current_val + random_shift
                new_val = max(min_val, min(max_val, new_val))
                
                param.setVal(new_val)
                param.setError(0.0)
        
        attempt += 1

    if not good_fit:
        print(f"WARNING: Fit failed to converge for {savename.stem} after {max_tries} attempts.")

    if attempt >= max_tries:
        print(f"WARNING: Fit failed to converge for {savename.stem} after {max_tries} attempts.")
        if get_params:
            return [None, None, None, None, None] if not signal_only else [None, None, None]
        return [None, None, None, None] if not signal_only else [None, None]


    check_parameter_limits(savename, fit_result)
    out = [sig_coeff.getVal(), sig_coeff.getError()]
    params_output = {
        'cb_mean'   : cb_mean.getVal(),
        'cb_sigma'  : cb_sigma.getVal(),
        'cb_alpha' : cb_alpha.getVal(),
        'cb_n'     : cb_n.getVal(),
    }

    if not signal_only:
        out.extend([bkg_coeff.getVal(), bkg_coeff.getError()])
        params_output.update({
            # **{i.GetName() : i for i in poly_pars},
            **{i.GetName() : i for i in bern_pars},
        })

    out += [params_output] if get_params else []

    frame = mass.frame(ROOT.RooFit.Title(str(savename.stem)))
    data.plotOn(frame)

    if signal_only:
        fit_model.plotOn(frame, ROOT.RooFit.Name(fit_model.GetName()), ROOT.RooFit.LineColor(38), ROOT.RooFit.Normalization(sig_coeff.getVal(),ROOT.RooAbsReal.NumEvent))
        h_pull = frame.pullHist()

    else:
        fit_model.plotOn(frame, ROOT.RooFit.Name(fit_model.GetName()), ROOT.RooFit.LineColor(38), ROOT.RooFit.Normalization(sig_coeff.getVal()+bkg_coeff.getVal(),ROOT.RooAbsReal.NumEvent))
        h_pull = frame.pullHist()
        fit_model.plotOn(frame, ROOT.RooFit.Components('sig_pdf'), ROOT.RooFit.LineColor(32), ROOT.RooFit.LineStyle(ROOT.kDashed))
        fit_model.plotOn(frame, ROOT.RooFit.Components('bkg_pdf'), ROOT.RooFit.LineColor(46), ROOT.RooFit.LineStyle(ROOT.kDashed))

    # fit_model.paramOn(frame, ROOT.RooFit.Layout(0.12, 0.3, 0.88), ROOT.RooFit.Format('NE', ROOT.RooFit.FixedPrecision(3)))

    frame_pull = mass.frame(ROOT.RooFit.Title(' '))
    frame_pull.addPlotable(h_pull, 'P')

    n_params = len(fit_result.floatParsFinal()) if fit_result else 0
    chi2 = frame.chiSquare('fit_model', 'h_data', n_params) if not signal_only else frame.chiSquare('fit_model', 'h_data', n_params)
    # chi2 = frame.chiSquare()
    chi2_text = ROOT.TLatex(0.7, 0.8, '#chi^{{2}}/ndf = {}'.format(round(chi2,1)))
    chi2_text.SetTextSize(0.05)
    chi2_text.SetNDC(ROOT.kTRUE)

    text = ROOT.TLatex(0.7, 0.7, f'N_{{J/#psi}} = {round(sig_coeff.getVal())} #pm {round(sig_coeff.getError(),1)}')
    text.SetTextSize(0.05)
    text.SetNDC(ROOT.kTRUE)

    c = ROOT.TCanvas('c', ' ', 800, 600)
    pad1 = ROOT.TPad('pad1', 'pad1', 0, 0.3, 1, 1.0)
    pad1.SetBottomMargin(0.02)
    pad1.SetGridx()
    pad1.Draw()
    c.cd()
    pad2 = ROOT.TPad('pad2', 'pad2', 0, 0.05, 1, 0.3)
    pad2.SetTopMargin(0.02)
    pad2.SetBottomMargin(0.2)
    pad2.SetGridx()
    pad2.Draw()

    pad1.cd()
    frame.Draw()
    ax_y_main = frame.GetYaxis()
    ax_x_main = frame.GetXaxis()
    ax_x_main.SetLabelOffset(3.)

    chi2_text.Draw()
    text.Draw()

    pad2.cd()
    frame_pull.Draw()

    ax_y_pull = frame_pull.GetYaxis()
    ax_x_pull = frame_pull.GetXaxis()

    line = ROOT.TLine(ax_x_pull.GetXmin(), 0, ax_x_pull.GetXmax(), 0)
    line.SetLineStyle(7)
    line.Draw()

    ax_y_pull.SetTitle('#frac{y - y_{fit}}{#sigma_{y}}')
    ax_y_pull.SetTitleOffset(.35)
    ax_y_pull.SetNdivisions(8)

    ax_y_pull.SetTitleSize(2.8*ax_y_main.GetTitleSize())
    ax_y_pull.SetLabelSize(2.8*ax_y_main.GetLabelSize())
    ax_x_pull.SetTitleSize(2.8*ax_x_main.GetTitleSize())
    ax_x_pull.SetLabelSize(2.8*ax_x_main.GetLabelSize())


    c.SaveAs(str(savename))
    c.Close()

    return out


class DotDict(dict):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        for key, value in self.items():
            if isinstance(value, dict):
                self[key] = DotDict(value)

    def __getattr__(self, key):
        return self.get(key, False)

    def __setattr__(self, key, value):
        self[key] = DotDict(value) if isinstance(value, dict) else value

    def __delattr__(self, key):
        if key in self: del self[key]
        else: raise AttributeError(f"No attribute named {key}")


def make_plotlist(cfg):
    plots = []
    out_path = Path('./eff_hists/test') if cfg.test else Path(cfg.output.output_dir)
    data_in_path = Path(cfg.inputs.data_dir).parent / 'test' if cfg.test else Path(cfg.inputs.data_dir)
    mc_in_path = Path(cfg.inputs.mc_dir).parent / 'test' if cfg.test else Path(cfg.inputs.mc_dir)
    for plot_name, plot_dict in cfg.plots.items():
        # For new format: triggers is a dict with data_trigger
        data_trigger = plot_dict.triggers['data_trigger'] if isinstance(plot_dict.triggers, dict) else plot_dict.triggers[0]
        for var in plot_dict.variables:
            plots.append(DotDict({
                'name': '_'.join([plot_name, f'{var}binned']),
                'data_trigger': data_trigger,
                'data_file': data_in_path / plot_dict.files.data,
                'mc_file': mc_in_path / plot_dict.files.mc if hasattr(plot_dict.files, 'mc') else None,
                'output_file': out_path / Path('_'.join(['eff', data_trigger, f'{var}binned'])).with_suffix('.png'),
                'mc_triggers': plot_dict.mc_triggers if hasattr(plot_dict, 'mc_triggers') else [],
                'var': var,
                **assign_hist_format(f'{var}binned'),
            }))
    return plots

def get_hists(cfg):
    # Only for data, since MC is handled per-path
    hists = []
    hist_paths = [(
        cfg.data_file, 
        f'hists/diel_m_{cfg.data_trigger}_num_{cfg.var}binned', 
        f'hists/diel_m_{cfg.data_trigger}_denom_{cfg.var}binned',
    )]
    for fname, num_path, denom_path in hist_paths:
        subhists = []
        f = ROOT.TFile(str(fname))
        num_h = f.Get(num_path)
        denom_h = f.Get(denom_path)
        nbins = num_h.GetYaxis().GetNbins()
        for i in range(1, nbins + 1):
            try:
                num_bin = num_h.ProjectionX('num_bin', i, i + 1)
                # print(f'bin {i} num: ',num_bin.GetEntries())
            except AttributeError:
                continue
            try:
                denom_bin = denom_h.ProjectionX('denom_bin', i, i + 1)
                # print(f'bin {i} deonom: ',denom_bin.GetEntries())
            except AttributeError:
                continue
            subhists.append((copy.deepcopy(num_bin), copy.deepcopy(denom_bin)))
        hists.append(subhists)
    return hists

def get_hists_2d(cfg):
    """Get 2D histograms for pt1pt2 binning"""
    hists = []
    hist_paths = [(
        cfg.data_file, 
        f'hists/diel_m_{cfg.data_trigger}_num_{cfg.var}binned', 
        f'hists/diel_m_{cfg.data_trigger}_denom_{cfg.var}binned',
    )]
    for fname, num_path, denom_path in hist_paths:
        subhists_2d = []
        f = ROOT.TFile(str(fname))
        num_h = f.Get(num_path)
        denom_h = f.Get(denom_path)
        
        # For 2D: iterate over Y bins (pt1), then X bins (pt2) for each pt1
        nbins_y = num_h.GetYaxis().GetNbins()  # pt1 bins
        nbins_z = num_h.GetZaxis().GetNbins()  # pt2 bins
        
        for iy in range(1, nbins_y + 1):
            pt1_bin_hists = []
            for iz in range(1, nbins_z + 1):
                try:
                    # Project to get mass distribution for this (pt1, pt2) bin
                    num_h.GetYaxis().SetRange(iy, iy)
                    num_h.GetZaxis().SetRange(iz, iz)
                    num_bin = num_h.Project3D("x")
                    num_bin.SetName(f'num_bin_pt1_{iy}_pt2_{iz}')
                except AttributeError:
                    continue
                try:
                    denom_h.GetYaxis().SetRange(iy, iy)
                    denom_h.GetZaxis().SetRange(iz, iz)
                    denom_bin = denom_h.Project3D("x")
                    denom_bin.SetName(f'denom_bin_pt1_{iy}_pt2_{iz}')
                except AttributeError:
                    continue
                pt1_bin_hists.append((copy.deepcopy(num_bin), copy.deepcopy(denom_bin)))
            subhists_2d.append(pt1_bin_hists)
        hists.append(subhists_2d)
    return hists

def make_eff_plot_dict(cfg):
    plot_list = make_plotlist(cfg)
    eff_plot_list = []
    for plot_cfg in plot_list:
        print(f'processing hists for {plot_cfg.name}')
        print("DEBUG: plot_cfg = ", plot_cfg)
        
        is_2d = plot_cfg.get('is_2d', False)
        
        eff_dict = {
            'name': plot_cfg.name,
            'data_trigger': plot_cfg.data_trigger,
            'output_file': plot_cfg.output_file,
            'data_num_yields': [] if not is_2d else [[]],
            'data_denom_yields': [] if not is_2d else [[]],
            'mc_eff_mixture': [] if not is_2d else [[]],
            'mc_triggers': plot_cfg.mc_triggers,
            'mc_name': plot_cfg.mc_file.name if plot_cfg.mc_file else None,
            'is_2d': is_2d,
        }
        
        if is_2d:
            eff_dict['bins_x'] = plot_cfg.bins_x
            eff_dict['bins_y'] = plot_cfg.bins_y
            eff_dict['xlabel'] = plot_cfg.xlabel
            eff_dict['ylabel'] = plot_cfg.ylabel
        else:
            eff_dict['bins'] = plot_cfg.bins
            eff_dict['xlabel'] = plot_cfg.xlabel
        
        # Data
        fit_output_file = plot_cfg.output_file.parent / 'fits' / plot_cfg.output_file.name
        
        if is_2d:
            # Handle 2D case
            data_hists_2d = get_hists_2d(plot_cfg)
            eff_dict['data_num_yields'] = []
            eff_dict['data_denom_yields'] = []
            
            for ipt1, pt1_bin_hists in enumerate(data_hists_2d[0]):
                num_yields_pt1 = []
                denom_yields_pt1 = []
                for ipt2, (data_num_hist, data_denom_hist) in enumerate(pt1_bin_hists):
                    num_results = do_fit(
                        data_num_hist,
                        savename=fit_output_file.with_stem(f'data_num_fit_{plot_cfg.name}_pt1bin{ipt1}_pt2bin{ipt2}'),
                        printlevel=cfg.printlevel,
                    )
                    n_num_data, n_num_data_err = (num_results[0], num_results[1]) if num_results[0] is not None else (None, None)

                    denom_results = do_fit(
                        data_denom_hist,
                        savename=fit_output_file.with_stem(f'data_denom_fit_{plot_cfg.name}_pt1bin{ipt1}_pt2bin{ipt2}'),
                        printlevel=cfg.printlevel,
                    )
                    n_denom_data, n_denom_data_err = (denom_results[0], denom_results[1]) if denom_results[0] is not None else (None, None)

                    num_yields_pt1.append((n_num_data, n_num_data_err))
                    denom_yields_pt1.append((n_denom_data, n_denom_data_err))
                
                eff_dict['data_num_yields'].append(num_yields_pt1)
                eff_dict['data_denom_yields'].append(denom_yields_pt1)
        else:
            # Handle 1D case (original code)
            data_hists = get_hists(plot_cfg)
            for i, (data_num_hist, data_denom_hist) in enumerate(data_hists[0]):
                num_results = do_fit(
                    data_num_hist,
                    savename=fit_output_file.with_stem(f'data_num_fit_{plot_cfg.name}_bin{i}'),
                    printlevel=cfg.printlevel,
                )
                n_num_data, n_num_data_err = (num_results[0], num_results[1]) if num_results[0] is not None else (None, None)

                denom_results = do_fit(
                    data_denom_hist,
                    savename=fit_output_file.with_stem(f'data_denom_fit_{plot_cfg.name}_bin{i}'),
                    printlevel=cfg.printlevel,
                )
                n_denom_data, n_denom_data_err = (denom_results[0], denom_results[1]) if denom_results[0] is not None else (None, None)

                eff_dict['data_num_yields'].append((n_num_data, n_num_data_err))
                eff_dict['data_denom_yields'].append((n_denom_data, n_denom_data_err))

        if plot_cfg.mc_file and plot_cfg.mc_triggers:
            mc_file = ROOT.TFile(str(plot_cfg.mc_file))
            
            if is_2d:
                # Handle 2D MC
                n_bins_x = len(eff_dict['bins_x']) - 1
                n_bins_y = len(eff_dict['bins_y']) - 1
                mc_path_effs = []
                
                for mc_trig in plot_cfg.mc_triggers:
                    path = mc_trig['path'].strip('"')
                    weight = mc_trig['weight']
                    num_hist_name = f'hists/diel_m_{path}_num_{plot_cfg.var}binned'
                    denom_hist_name = f'hists/diel_m_{path}_denom_{plot_cfg.var}binned'
                    num_hist = mc_file.Get(num_hist_name)
                    denom_hist = mc_file.Get(denom_hist_name)
                    
                    effs_2d = []
                    for ipt1 in range(1, n_bins_y + 1):
                        effs_pt1 = []
                        for ipt2 in range(1, n_bins_x + 1):
                            num_hist.GetYaxis().SetRange(ipt1, ipt1)
                            num_hist.GetZaxis().SetRange(ipt2, ipt2)
                            num_proj = num_hist.Project3D("x")
                            
                            denom_hist.GetYaxis().SetRange(ipt1, ipt1)
                            denom_hist.GetZaxis().SetRange(ipt2, ipt2)
                            denom_proj = denom_hist.Project3D("x")
                            
                            if denom_proj.GetSumOfWeights() > 0:
                                nNum, nNum_err = get_event_count(num_proj)
                                nDenom, nDenom_err = get_event_count(denom_proj)
                                eff = ufloat(nNum, nNum_err) / ufloat(nDenom, nDenom_err)
                            else:
                                eff = ufloat(0, 0)
                            effs_pt1.append(eff * weight)
                        effs_2d.append(effs_pt1)
                    mc_path_effs.append(effs_2d)

                mc_eff_mixture = []
                for ipt1 in range(n_bins_y):
                    mixture_pt1 = []
                    for ipt2 in range(n_bins_x):
                        eff_sum = sum(effs[ipt1][ipt2] for effs in mc_path_effs)
                        mixture_pt1.append(eff_sum)
                    mc_eff_mixture.append(mixture_pt1)
                eff_dict['mc_eff_mixture'] = mc_eff_mixture
            else:
                # Handle 1D MC (original code)
                n_bins = len(eff_dict['bins']) - 1
                mc_path_effs = []
                for mc_trig in plot_cfg.mc_triggers:
                    path = mc_trig['path'].strip('"')
                    weight = mc_trig['weight']
                    num_hist_name = f'hists/diel_m_{path}_num_{plot_cfg.var}binned'
                    denom_hist_name = f'hists/diel_m_{path}_denom_{plot_cfg.var}binned'
                    num_hist = mc_file.Get(num_hist_name)
                    denom_hist = mc_file.Get(denom_hist_name)
                    effs = []
                    for ibin in range(1, n_bins + 1):
                        num_proj = num_hist.ProjectionX('num_bin', ibin, ibin + 1)
                        denom_proj = denom_hist.ProjectionX('denom_bin', ibin, ibin + 1)
                        
                        if denom_proj.GetSumOfWeights() > 0:
                            nNum, nNum_err = get_event_count(num_proj)
                            nDenom, nDenom_err = get_event_count(denom_proj)
                            eff = ufloat(nNum, nNum_err) / ufloat(nDenom, nDenom_err)
                        else:
                            eff = ufloat(0, 0)
                        print("DEBUG: eff = ", eff, "; weight = ", weight)
                        effs.append(eff * weight)
                        print(f"mc_trig {path}: {weight}, eff: {eff}, eff*weight: {eff * weight}, num: {nNum}, denom: {nDenom}")
                    mc_path_effs.append(effs)

                mc_eff_mixture = []
                for ibin in range(n_bins):
                    eff_sum = sum(effs[ibin] for effs in mc_path_effs)
                    mc_eff_mixture.append(eff_sum)
                eff_dict['mc_eff_mixture'] = mc_eff_mixture
        eff_plot_list.append(eff_dict)

    with open(Path(cfg.output.output_dir) / 'eff_plot_entries.pkl', 'wb') as f:
        pkl.dump(eff_plot_list, f)
    
    return eff_plot_list


def plot_efficiencies(eff_dicts, test=False):
    for d in eff_dicts:
        if test and ('test' not in d['name']):
            continue
        elif not test and ('test' in d['name']):
            continue

        if d.get('is_2d', False):
            # Handle 2D case
            plot_efficiencies_2d(d)
        else:
            # Handle 1D case (original code)
            bins = d['bins']
            h_num_data = ROOT.TH1F('h_num_data', 'h_num_data', len(bins)-1, bins)
            h_denom_data = ROOT.TH1F('h_denom_data', 'h_denom_data', len(bins)-1, bins)
            h_num_mc = ROOT.TH1F('h_num_mc', 'h_num_mc', len(bins)-1, bins)
            h_denom_mc = ROOT.TH1F('h_denom_mc', 'h_denom_mc', len(bins)-1, bins)

            # Fill data as before
            for ibin, (data_num, data_denom) in enumerate(zip(d['data_num_yields'], d['data_denom_yields'])):
                if data_num[0] is not None and data_denom[0] is not None:
                    h_num_data.SetBinContent(ibin+1, data_num[0])
                    h_num_data.SetBinError(ibin+1, data_num[1])
                    h_denom_data.SetBinContent(ibin+1, data_denom[0])
                    h_denom_data.SetBinError(ibin+1, data_denom[1])


            # MC: use mixture if present
            if 'mc_eff_mixture' in d and d['mc_eff_mixture']:
                for ibin, eff in enumerate(d['mc_eff_mixture']):
                    # For TEfficiency, set numerator = eff * denom, error propagation is handled by ufloat
                    denom = 1.0  # or set to 1 for all bins, since TEfficiency cares about ratio
                    num = eff.n
                    err = eff.s
                    h_num_mc.SetBinContent(ibin+1, num)
                    h_num_mc.SetBinError(ibin+1, err)
                    h_denom_mc.SetBinContent(ibin+1, denom)
                    h_denom_mc.SetBinError(ibin+1, 0)
            else:
                for ibin, (mc_num, mc_denom) in enumerate(zip(d['mc_num_yields'], d['mc_denom_yields'])):
                    h_num_mc.SetBinContent(ibin+1, mc_num[0])
                    h_num_mc.SetBinError(ibin+1, mc_num[1])
                    h_denom_mc.SetBinContent(ibin+1, mc_denom[0])
                    h_denom_mc.SetBinError(ibin+1, mc_denom[1])

            eff_data = ROOT.TEfficiency(h_num_data, h_denom_data)
            eff_data.SetStatisticOption(ROOT.TEfficiency.kBBayesian)
            eff_mc = ROOT.TEfficiency(h_num_mc, h_denom_mc)
            eff_mc.SetStatisticOption(ROOT.TEfficiency.kFNormal)

            eff_plot = EfficiencyPlot(init_params={
                'title_string': ' ;' + d["xlabel"] + ';Efficiency',
                'xrange': (bins[0], bins[-1]),
                'yrange': (0., 1.1),
                'leg_scale': .65,
                'leg_header': '',
                'x_title': d["xlabel"],
            })

            eff_plot.plotEfficiencies(
                eff_data,
                eff_mc,
                ratio=True,
                h1_title='ParkingDoubleMuonLowMass Data',
                h2_title='Prompt J/#psi MC' if 'prompt' in d['mc_name'].lower() else 'B^{+} #rightarrow J/#psi K^{+} MC',
                save=str(d['output_file']),
                addIntegral=False
            )


def plot_efficiencies_2d(d):
    """Plot 2D efficiency maps"""
    bins_x = d['bins_x']
    bins_y = d['bins_y']
    
    # Create 2D histograms for data
    h_eff_data = ROOT.TH2F('h_eff_data', 'Data Efficiency;' + d['xlabel'] + ';' + d['ylabel'], 
                           len(bins_x)-1, bins_x, len(bins_y)-1, bins_y)
    h_eff_mc = ROOT.TH2F('h_eff_mc', 'MC Efficiency;' + d['xlabel'] + ';' + d['ylabel'], 
                         len(bins_x)-1, bins_x, len(bins_y)-1, bins_y)
    h_sf = ROOT.TH2F('h_sf', 'Scale Factor;' + d['xlabel'] + ';' + d['ylabel'], 
                     len(bins_x)-1, bins_x, len(bins_y)-1, bins_y)
    
    # Fill data efficiency
    # ipt1 corresponds to Y-axis (leading pT) in the 3D histogram
    # ipt2 corresponds to Z-axis (subleading pT) in the 3D histogram
    # For 2D plot: X-axis = leading (pt1), Y-axis = subleading (pt2)
    for ipt1, (num_yields_pt1, denom_yields_pt1) in enumerate(zip(d['data_num_yields'], d['data_denom_yields'])):
        for ipt2, (data_num, data_denom) in enumerate(zip(num_yields_pt1, denom_yields_pt1)):
            if data_num[0] is not None and data_denom[0] is not None and data_denom[0] > 0:
                eff = data_num[0] / data_denom[0]
                # Error propagation for efficiency
                err = eff * np.sqrt((data_num[1]/data_num[0])**2 + (data_denom[1]/data_denom[0])**2) if data_num[0] > 0 else 0
                # SetBinContent(x_bin, y_bin, value): x = leading (ipt1), y = subleading (ipt2)
                h_eff_data.SetBinContent(ipt1+1, ipt2+1, eff)
                h_eff_data.SetBinError(ipt1+1, ipt2+1, err)
    
    # Fill MC efficiency
    if 'mc_eff_mixture' in d and d['mc_eff_mixture']:
        for ipt1, effs_pt1 in enumerate(d['mc_eff_mixture']):
            for ipt2, eff in enumerate(effs_pt1):
                h_eff_mc.SetBinContent(ipt1+1, ipt2+1, eff.n)
                h_eff_mc.SetBinError(ipt1+1, ipt2+1, eff.s)
    
    # Calculate scale factors
    for ipt1 in range(1, len(bins_x)):
        for ipt2 in range(1, len(bins_y)):
            eff_data = h_eff_data.GetBinContent(ipt1, ipt2)
            eff_mc = h_eff_mc.GetBinContent(ipt1, ipt2)
            err_data = h_eff_data.GetBinError(ipt1, ipt2)
            err_mc = h_eff_mc.GetBinError(ipt1, ipt2)
            
            if eff_mc > 0:
                sf = eff_data / eff_mc
                # Error propagation for ratio
                err_sf = sf * np.sqrt((err_data/eff_data)**2 + (err_mc/eff_mc)**2) if eff_data > 0 else 0
                h_sf.SetBinContent(ipt1, ipt2, sf)
                h_sf.SetBinError(ipt1, ipt2, err_sf)
    
    # Create output plots
    ROOT.gStyle.SetOptStat(0)
    ROOT.gStyle.SetPaintTextFormat("0.3f")
    
    for hist, suffix, title in [(h_eff_data, 'data_eff', 'Data Efficiency'), 
                                 (h_eff_mc, 'mc_eff', 'MC Efficiency'), 
                                 (h_sf, 'sf', 'Scale Factor')]:
        c = ROOT.TCanvas('c_' + suffix, title, 800, 600)
        c.SetRightMargin(0.15)
        hist.SetTitle(title)
        hist.SetMinimum(0)
        hist.SetMaximum(1.0 if 'sf' not in suffix else 1.5)
        hist.Draw('COLZ TEXT')
        
        output_file = d['output_file'].with_stem(d['output_file'].stem + '_' + suffix)
        c.SaveAs(str(output_file))
        c.Close()
    
    # Now create sliced 1D plots
    plot_efficiencies_2d_slices(d, h_eff_data, h_eff_mc, bins_x, bins_y)


def plot_efficiencies_2d_slices(d, h_eff_data, h_eff_mc, bins_x, bins_y):
    """Create 1D sliced efficiency plots from 2D histogram"""
    
    # Slice by leading pT: for each leading pT bin, plot efficiency vs subleading pT
    for ipt1 in range(1, len(bins_x)):
        h_data_slice = ROOT.TH1F(f'h_data_leadpt{ipt1}', 
                                 f'Data Eff, {bins_x[ipt1-1]:.0f} < p_{{T,lead}} < {bins_x[ipt1]:.0f} GeV',
                                 len(bins_y)-1, bins_y)
        h_mc_slice = ROOT.TH1F(f'h_mc_leadpt{ipt1}', 
                               f'MC Eff, {bins_x[ipt1-1]:.0f} < p_{{T,lead}} < {bins_x[ipt1]:.0f} GeV',
                               len(bins_y)-1, bins_y)
        
        for ipt2 in range(1, len(bins_y)):
            # Extract data and MC values for this leading pT bin
            eff_data = h_eff_data.GetBinContent(ipt1, ipt2)
            err_data = h_eff_data.GetBinError(ipt1, ipt2)
            eff_mc = h_eff_mc.GetBinContent(ipt1, ipt2)
            err_mc = h_eff_mc.GetBinError(ipt1, ipt2)
            
            h_data_slice.SetBinContent(ipt2, eff_data)
            h_data_slice.SetBinError(ipt2, err_data)
            h_mc_slice.SetBinContent(ipt2, eff_mc)
            h_mc_slice.SetBinError(ipt2, err_mc)
        
        # Plot this slice
        c = ROOT.TCanvas(f'c_slice_leadpt{ipt1}', f'Slice Leading pT bin {ipt1}', 800, 600)
        c.SetGridy()
        
        h_data_slice.SetLineColor(ROOT.kBlack)
        h_data_slice.SetMarkerStyle(20)
        h_data_slice.SetMarkerColor(ROOT.kBlack)
        h_data_slice.GetYaxis().SetTitle('Efficiency')
        h_data_slice.GetXaxis().SetTitle(d['ylabel'])
        h_data_slice.SetMinimum(0)
        h_data_slice.SetMaximum(1.1)
        h_data_slice.Draw('E1')
        
        h_mc_slice.SetLineColor(ROOT.kRed)
        h_mc_slice.SetMarkerStyle(24)
        h_mc_slice.SetMarkerColor(ROOT.kRed)
        h_mc_slice.Draw('E1 SAME')
        
        leg = ROOT.TLegend(0.6, 0.15, 0.88, 0.35)
        leg.AddEntry(h_data_slice, 'Data', 'lep')
        leg.AddEntry(h_mc_slice, 'MC', 'lep')
        leg.Draw()
        
        output_file = d['output_file'].with_stem(d['output_file'].stem + f'_slice_leadpt_bin{ipt1}')
        c.SaveAs(str(output_file))
        c.Close()
    
    # Slice by subleading pT: for each subleading pT bin, plot efficiency vs leading pT
    for ipt2 in range(1, len(bins_y)):
        h_data_slice = ROOT.TH1F(f'h_data_subpt{ipt2}', 
                                 f'Data Eff, {bins_y[ipt2-1]:.0f} < p_{{T,sublead}} < {bins_y[ipt2]:.0f} GeV',
                                 len(bins_x)-1, bins_x)
        h_mc_slice = ROOT.TH1F(f'h_mc_subpt{ipt2}', 
                               f'MC Eff, {bins_y[ipt2-1]:.0f} < p_{{T,sublead}} < {bins_y[ipt2]:.0f} GeV',
                               len(bins_x)-1, bins_x)
        
        for ipt1 in range(1, len(bins_x)):
            # Extract data and MC values for this subleading pT bin
            eff_data = h_eff_data.GetBinContent(ipt1, ipt2)
            err_data = h_eff_data.GetBinError(ipt1, ipt2)
            eff_mc = h_eff_mc.GetBinContent(ipt1, ipt2)
            err_mc = h_eff_mc.GetBinError(ipt1, ipt2)
            
            h_data_slice.SetBinContent(ipt1, eff_data)
            h_data_slice.SetBinError(ipt1, err_data)
            h_mc_slice.SetBinContent(ipt1, eff_mc)
            h_mc_slice.SetBinError(ipt1, err_mc)
        
        # Plot this slice
        c = ROOT.TCanvas(f'c_slice_subpt{ipt2}', f'Slice Subleading pT bin {ipt2}', 800, 600)
        c.SetGridy()
        
        h_data_slice.SetLineColor(ROOT.kBlack)
        h_data_slice.SetMarkerStyle(20)
        h_data_slice.SetMarkerColor(ROOT.kBlack)
        h_data_slice.GetYaxis().SetTitle('Efficiency')
        h_data_slice.GetXaxis().SetTitle(d['xlabel'])
        h_data_slice.SetMinimum(0)
        h_data_slice.SetMaximum(1.1)
        h_data_slice.Draw('E1')
        
        h_mc_slice.SetLineColor(ROOT.kRed)
        h_mc_slice.SetMarkerStyle(24)
        h_mc_slice.SetMarkerColor(ROOT.kRed)
        h_mc_slice.Draw('E1 SAME')
        
        leg = ROOT.TLegend(0.6, 0.15, 0.88, 0.35)
        leg.AddEntry(h_data_slice, 'Data', 'lep')
        leg.AddEntry(h_mc_slice, 'MC', 'lep')
        leg.Draw()
        
        output_file = d['output_file'].with_stem(d['output_file'].stem + f'_slice_subpt_bin{ipt2}')
        c.SaveAs(str(output_file))
        c.Close()


def make_sf_json(cfg, eff_dicts):
    outputs = {}
    for eff_dict in eff_dicts:
        matches = list(re.finditer(r'_([^_]+)binned(?=_|$)', eff_dict['name']))
        if not matches:
            raise ValueError(f"Could not find binvar in {eff_dict['name']}")

        binvar = matches[-1].group(1)
        
        if eff_dict.get('is_2d', False):
            # Handle 2D case
            data_num_2d = [[ufloat(val, unc) if val is not None else None for val, unc in row] 
                          for row in eff_dict['data_num_yields']]
            data_denom_2d = [[ufloat(val, unc) if val is not None else None for val, unc in row] 
                            for row in eff_dict['data_denom_yields']]
            
            # Calculate data ratios
            data_ratio_2d = []
            for num_row, denom_row in zip(data_num_2d, data_denom_2d):
                ratio_row = [(n / d) if (n is not None and d is not None and d.n != 0) else None 
                            for n, d in zip(num_row, denom_row)]
                data_ratio_2d.append(ratio_row)
            
            mc_ratio_2d = eff_dict['mc_eff_mixture']
            
            # Calculate SFs
            sfs_2d = []
            for data_row, mc_row in zip(data_ratio_2d, mc_ratio_2d):
                sf_row = [(d / m) if (d is not None and m.n != 0) else None 
                         for d, m in zip(data_row, mc_row)]
                sfs_2d.append(sf_row)
            
            # Format output
            sfs_out = [[(round(r.n, 3), round(r.s, 3)) if r is not None else None for r in row] 
                      for row in sfs_2d]
            data_ratio_out = [[(round(r.n, 3), round(r.s, 3)) if r is not None else None for r in row] 
                             for row in data_ratio_2d]
            mc_ratio_out = [[(round(r.n, 3), round(r.s, 3)) if r is not None else None for r in row] 
                           for row in mc_ratio_2d]
            
            output = {
                'binvar'    : binvar,
                'bins_x'    : list(eff_dict['bins_x']),
                'bins_y'    : list(eff_dict['bins_y']),
                'sfs'       : sfs_out,
                'data_effs' : data_ratio_out,
                'mc_effs'   : mc_ratio_out,
                'is_2d'     : True,
                'sfs_ufloat': sfs_2d,  # Keep ufloat for ROOT file
                'data_effs_ufloat': data_ratio_2d,
                'mc_effs_ufloat': mc_ratio_2d,
            }
        else:
            # Handle 1D case (original code)
            data_num = [ufloat(val, unc) if val is not None else None for val, unc in eff_dict['data_num_yields']]
            data_denom = [ufloat(val, unc) if val is not None else None for val, unc in eff_dict['data_denom_yields']]

            # Propagate None through ratios
            data_ratio = [(n / d) if (n is not None and d is not None and d.n != 0) else None for n, d in zip(data_num, data_denom)]
            mc_ratio = eff_dict['mc_eff_mixture']
            sfs = [(d / m) if (d is not None and m.n != 0) else None for d, m in zip(data_ratio, mc_ratio)]

            # Format final output, preserving None
            sfs_out = [ (round(r.n, 3), round(r.s, 3)) if r is not None else None for r in sfs]
            data_ratio_out = [ (round(r.n, 3), round(r.s, 3)) if r is not None else None for r in data_ratio]
            mc_ratio_out = [ (round(r.n, 3), round(r.s, 3)) if r is not None else None for r in mc_ratio]

            output = {
                'binvar'    : binvar,
                'bins'      : list(eff_dict['bins']),
                'sfs'       : sfs_out,
                'data_effs' : data_ratio_out,
                'mc_effs'   : mc_ratio_out,
                'is_2d'     : False,
                'sfs_ufloat': sfs,  # Keep ufloat for ROOT file
                'data_effs_ufloat': data_ratio,
                'mc_effs_ufloat': mc_ratio,
            }

        outputs[eff_dict['name']] = outputs.get(eff_dict['name'], copy.deepcopy(output))

    out_path = Path('eff_hists') / 'test' if cfg.test else Path(cfg.output.output_dir)
    
    # Save JSON
    with open(out_path / 'sf_jsons' / 'trigger_sfs.json', 'w') as outfile:
        json.dump({k: {key: val for key, val in v.items() if 'ufloat' not in key} for k, v in outputs.items()}, outfile, indent=4)
    
    # Save ROOT file
    make_sf_root_file(cfg, outputs)


def make_sf_root_file(cfg, outputs):
    """Create ROOT file with SF histograms"""
    out_path = Path('eff_hists') / 'test' if cfg.test else Path(cfg.output.output_dir)
    root_file = ROOT.TFile(str(out_path / 'sf_jsons' / 'trigger_sfs.root'), 'RECREATE')
    
    for name, output in outputs.items():
        if output['is_2d']:
            # Create 2D histograms
            bins_x = np.array(output['bins_x'], dtype=np.double)
            bins_y = np.array(output['bins_y'], dtype=np.double)
            
            h_sf = ROOT.TH2F(f'sf_{name}', f'Scale Factor {name};Leading p_{{T}} [GeV];Subleading p_{{T}} [GeV]',
                            len(bins_x)-1, bins_x, len(bins_y)-1, bins_y)
            h_data_eff = ROOT.TH2F(f'data_eff_{name}', f'Data Efficiency {name};Leading p_{{T}} [GeV];Subleading p_{{T}} [GeV]',
                                  len(bins_x)-1, bins_x, len(bins_y)-1, bins_y)
            h_mc_eff = ROOT.TH2F(f'mc_eff_{name}', f'MC Efficiency {name};Leading p_{{T}} [GeV];Subleading p_{{T}} [GeV]',
                                len(bins_x)-1, bins_x, len(bins_y)-1, bins_y)
            
            for ipt1, (sf_row, data_row, mc_row) in enumerate(zip(output['sfs_ufloat'], output['data_effs_ufloat'], output['mc_effs_ufloat'])):
                for ipt2, (sf, data_eff, mc_eff) in enumerate(zip(sf_row, data_row, mc_row)):
                    if sf is not None:
                        h_sf.SetBinContent(ipt1+1, ipt2+1, sf.n)
                        h_sf.SetBinError(ipt1+1, ipt2+1, sf.s)
                    if data_eff is not None:
                        h_data_eff.SetBinContent(ipt1+1, ipt2+1, data_eff.n)
                        h_data_eff.SetBinError(ipt1+1, ipt2+1, data_eff.s)
                    if mc_eff is not None:
                        h_mc_eff.SetBinContent(ipt1+1, ipt2+1, mc_eff.n)
                        h_mc_eff.SetBinError(ipt1+1, ipt2+1, mc_eff.s)
            
            h_sf.Write()
            h_data_eff.Write()
            h_mc_eff.Write()
        else:
            # Create 1D histograms
            bins = np.array(output['bins'], dtype=np.double)
            
            h_sf = ROOT.TH1F(f'sf_{name}', f'Scale Factor {name};p_{{T}} [GeV];SF',
                            len(bins)-1, bins)
            h_data_eff = ROOT.TH1F(f'data_eff_{name}', f'Data Efficiency {name};p_{{T}} [GeV];Efficiency',
                                  len(bins)-1, bins)
            h_mc_eff = ROOT.TH1F(f'mc_eff_{name}', f'MC Efficiency {name};p_{{T}} [GeV];Efficiency',
                                len(bins)-1, bins)
            
            for ibin, (sf, data_eff, mc_eff) in enumerate(zip(output['sfs_ufloat'], output['data_effs_ufloat'], output['mc_effs_ufloat'])):
                if sf is not None:
                    h_sf.SetBinContent(ibin+1, sf.n)
                    h_sf.SetBinError(ibin+1, sf.s)
                if data_eff is not None:
                    h_data_eff.SetBinContent(ibin+1, data_eff.n)
                    h_data_eff.SetBinError(ibin+1, data_eff.s)
                if mc_eff is not None:
                    h_mc_eff.SetBinContent(ibin+1, mc_eff.n)
                    h_mc_eff.SetBinError(ibin+1, mc_eff.s)
            
            h_sf.Write()
            h_data_eff.Write()
            h_mc_eff.Write()
    
    root_file.Close()
    print(f"ROOT file saved to {out_path / 'sf_jsons' / 'trigger_sfs.root'}")

def main(cfg):
    cfg = DotDict(cfg)

    if cfg.file is None:
        eff_dicts = make_eff_plot_dict(cfg)
    elif Path(cfg.file).is_file():
        with open(cfg.file,'rb') as f:
            eff_dicts = pkl.load(f)
    else:
        raise ValueError(f'Efficiency dict file {str(cfg.file)} is not readable')

    plot_efficiencies(eff_dicts, test=cfg.test)
    make_sf_json(cfg, eff_dicts)


if __name__=='__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-c', '--config', dest='config', type=str, default='eff_plot_cfg_v2.yml', help='plot configuration file (.yml)')
    parser.add_argument('-v', '--verbose', dest='verbose', action='store_true', help='printouts to stdout')
    parser.add_argument('-t', '--test', dest='test', action='store_true', help='only run test samples')
    parser.add_argument('-f', '--file', dest='file', default=None, help='make plots from pkl file')
    args = parser.parse_args()

    with open(args.config, 'r') as f:
        cfg = DotDict(yaml.safe_load(f))

    os.makedirs(Path(cfg.output.output_dir), exist_ok=True)
    os.makedirs(Path(cfg.output.output_dir) / 'fits', exist_ok=True)
    os.makedirs(Path(cfg.output.output_dir) / 'sf_jsons', exist_ok=True)
    cfg.test = args.test
    cfg.verbose = args.verbose
    cfg.printlevel = set_verbosity(args.verbose)
    cfg.file = args.file

    main(cfg)

