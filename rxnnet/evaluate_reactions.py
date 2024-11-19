#  This file is part of RxnNet.
#  Copyright © 2024 Technical University of Denmark

import os.path

import numpy as np
import pandas as pd
import wquantiles

pd.options.mode.chained_assignment = None


def weighted_median(target, target_weight):
    """
    Calculates the weighted median.

    :param target: target values
    :param target_weight: target weight values
    :return: weighted median (float)
    """
    mask = ~target.isna() & ~target_weight.isna()
    target, target_weight = target[mask], target_weight[mask]
    return wquantiles.median(target, target_weight)


def structuralsimilarity_filter(rn_eval) -> pd.DataFrame:
    """
    Calculate structural similarity weights from reaction network dataframe w/ reaction-based structural similarity distances.

    :param rn_eval: processed RN dataframe (pandas Dataframe)
    :return: processed RN dataframe with SSW weights (pandas Dataframe)
    """
    w_abs = rn_eval['d_rxn']
    w_abs *= w_abs.size / np.sum(w_abs)  # normalize
    rn_eval['weight'] *= np.max(w_abs) - w_abs  # relative weight
    df_rxn = rn_eval[rn_eval['weight'] != 0]  # remove reaction with zero weight
    return df_rxn


class EvaluateRN:
    """
    Class used to evaluate and obtain predictions from reaction files.
    """

    def __init__(self, target_id: str, preferred_method: str, reference_method: str, rn: str,
                 db=r'..\rxnnet\data\database-20.03.24\datasets\main.pkl.gz', mode='ssw+cf', fallback_method=None,
                 filter_cutoff=10, temperature=None, weight_stoic=False) -> None:
        """
            :param target_id: target compound identifier (str)
            :param preferred_method: alias of preferred calculation method (str)
            :param reference_method: alias of reference method (str)
            :param rn: filename of reaction network dataframe for target compound (str)
            :param db: filename of database containing at minimum the calculated and reference values of compounds in reaction network (str)
            :param mode: reaction weighing and selection mode to be used to attempt to improve error-cancellation [choose from: 'plain', 'ssw', 'cf' and 'ssw+cf' (default)] (str)
            :param fallback_method: alias of calculation method to be used in case a reaction cannot be composed using data purely from the preferred method (str)
            :param filter_cutoff: minimum number of reactions to allow for prediction if mode is not 'plain' (int)
            :param temperature: if input data is provided as dictionary of "temperature,value"-format a temperature may be specified. (dict)
            :param weight_stoic: enable a weighing mechanism which introduces a crude accounting of the random error in the calculated data (bool)
        """
        self.target_id = target_id
        self.preferred_method = preferred_method
        self.reference_method = reference_method
        self.db = pd.read_pickle(db)
        df_rn = pd.read_pickle(rn)
        df_rn['weight'] = np.ones(len(df_rn.index))  # set default weights (will be overwritten if ssw in mode).
        self.rn = df_rn
        self.mode = mode
        self.fallback_method = fallback_method
        self.filter_cutoff = filter_cutoff
        self.temperature = temperature
        self.rn_out = None
        self.weight_stoic = weight_stoic

    def get_prediction(self, reaction_ids: list[list, list], stoichiometry: list[list, list]):

        """
        Generate prediction for a single reaction, taking into account whether a hybrid method is being used.

        :param reaction_ids: ids of reactants and products (nested list)
        :param stoichiometry:  stoichiometry of reactants and products (nested list)
        :return: predicted value and calculation method used to make prediction (tuple)
        """

        target_id = self.target_id
        preferred_method = self.preferred_method
        db = self.db
        fallback_method = self.fallback_method
        reference_method = self.reference_method
        temperature = self.temperature

        def get_prediction_value(select_method: str):
            """
            Calculate prediction using a given method.

            :param select_method: alias of calculation method (str)
            :return: predicted value (float)
            """

            def get_values(compound_id):
                """
                Obtain calculated and experimental values for specific compound.
                :param compound_id: compound identifier
                :return: Calculated and experimental value for specific compound (tuple).
                """
                calc = db.loc[compound_id, select_method]
                exp = db.loc[compound_id, reference_method]
                if temperature:
                    if calc == calc:
                        calc = calc[temperature]
                    if exp == exp:
                        exp = exp[temperature]
                return calc, exp

            e_rxn_calc, e_rxn_ref, v0 = 0, 0, 1
            for side in reaction_ids:
                for i, c_id in enumerate(reaction_ids[side]):
                    is_target_compound = c_id == target_id
                    e_calc, e_ref = get_values(c_id)
                    if e_calc is None or (e_ref is None and not is_target_compound):
                        return np.nan, np.nan
                    v = stoichiometry[side][i]
                    e_rxn_calc += e_calc * v
                    if is_target_compound:
                        v0 = abs(v)
                    else:
                        e_rxn_ref += e_ref * v
            return (e_rxn_ref - e_rxn_calc) / v0, select_method

        if preferred_method is None:
            preferred_method = fallback_method

        if preferred_method != fallback_method:
            pred, out_method = get_prediction_value(preferred_method)
            if pd.isna(pred) and not pd.isna(fallback_method):
                pred, out_method = get_prediction_value(fallback_method)
        else:
            pred, out_method = get_prediction_value(fallback_method)
        return pred, out_method

    def get_evaluated_reactions(self, exclude_id=()):
        """
        For a compound, obtain predictions for every reaction in its reaction network.
        :param exclude_id: list of compounds which should be excluded from the reaction network (tuple)
        :return: dataframe of reactions with predictions (pandas Dataframe)
        """
        db = self.db
        rn = self.rn
        weight_stoic = self.weight_stoic

        if len(rn.index) == 0:
            return rn
        id_flattened = tuple(
            [[mpid for mpid_tuple in x.values() for mpid in mpid_tuple] for i, x in enumerate(rn['reaction_mpids'])])
        rn = rn[[set(c_ids).issubset(set(db.index.to_list()) - set(exclude_id)) for c_ids in id_flattened]]
        rn.loc[:, ['value', 'method']] = [self.get_prediction(x, rn['stoichiometry (atoms)'].iloc[i]) for i, x in
                                          enumerate(rn['reaction_mpids'])]
        if weight_stoic:
            stoic_weights = []
            for brxn in rn['BalancedReaction']:
                err, v0 = 0, 1
                for i, r in enumerate(brxn['reactants'].values()):
                    err += r * r
                    if i == 0:
                        v0 = r
                for p in brxn['products'].values():
                    err += p * p
                stoic_weights.append(1 / (np.sqrt(err) / v0))
            stoic_weights = np.array(stoic_weights)
            rn['weight'] *= stoic_weights / np.sum(stoic_weights)
        if len(rn.index) == 0:
            return rn
        rn.dropna(inplace=True)
        return pd.DataFrame(rn)

    def process_reaction_mode(self, rn_eval: pd.DataFrame, mode: str):
        """
        Depending on error-cancellation mode evaluate and re-evaluate RN dataframe accordingly to return RN dataframe
        containing optimal reactions.

        :param rn_eval: processed RN dataframe (pd.Dataframe)
        :param mode: error-cancellation mode to be used ('ssw', 'cf', 'ssw+cf')
        :return: processed RN dataframe (pd.Dataframe)
        """
        rn_eval_in = rn_eval.copy()
        cutoff = self.filter_cutoff

        if rn_eval.empty:  # empty or non-existing dataframe
            return None, None
        elif mode == 'plain':
            return rn_eval, mode
        elif mode == 'ssw':
            rn_eval = structuralsimilarity_filter(rn_eval)
            if len(rn_eval.index) <= cutoff:
                rn_eval, mode = self.process_reaction_mode(rn_eval_in, mode='plain')
                return rn_eval, mode
            else:
                return rn_eval, mode
        elif mode == 'cf':
            rn_eval = rn_eval[rn_eval['balanced_chemistry']]
            if len(rn_eval.index) <= cutoff:
                rn_eval, mode = self.process_reaction_mode(rn_eval_in, mode='plain')
                return rn_eval, mode
            else:
                return rn_eval, mode
        elif mode == 'ssw+cf':
            rn_eval = rn_eval[rn_eval['balanced_chemistry']]
            if len(rn_eval.index) <= cutoff + 1:  # ssw will remove an additional reaction
                rn_eval, mode = self.process_reaction_mode(rn_eval_in, mode='ssw')
                return rn_eval, mode
            else:
                rn_eval = structuralsimilarity_filter(rn_eval)
                return rn_eval, mode
        else:
            return None, None

    def rn_evaluate(self, save=False, verbose=False, save_path='processed_reactions', exclude_id=()):
        """
        Computes the final RN prediction along with uncertainty information,
        details about the reaction distribution/network and general information about the compound.

        :param verbose: whether to print table with detailed information about reactions and predictions (bool).
        :param save: whether to save dataframe with detailed information about reactions and predictions (bool).
        :param save_path: where to save dataframe with detailed information about reactions and predictions (str).
        :param exclude_id: list of compounds which should be excluded from the reaction network (tuple)

        :return: summary data pertaining to and including the final RN prediction(s) (dict)
        """

        target_id = self.target_id
        preferred_method = self.preferred_method
        fallback_method = self.fallback_method
        db = self.db
        mode = self.mode

        def wmean_wmed(df):
            """

            :param df:
            :return:
            """
            if len(df.index) > 0:
                wmean, wmed = np.ma.average(np.ma.masked_invalid(df['value']), weights=df['weight']), weighted_median(
                    df['value'], df['weight'])
            else:
                wmean, wmed = np.nan, np.nan
            return wmean, wmed

        data = []
        if target_id in exclude_id:
            exclude_list = list(exclude_id)
            exclude_list.remove(target_id)
            exclude_id = tuple(exclude_list)
        if mode not in ['plain', 'ssw', 'cf', 'ssw+cf']:
            print("Mode not recognized. Using mode 'ssw+cf'. Supported modes are: 'plain', 'ssw', 'cf', 'ssw+cf'.")
            mode = 'ssw+cf'
        rn_eval = self.get_evaluated_reactions(exclude_id)
        df_out, filter_type = self.process_reaction_mode(rn_eval, mode)
        if df_out is None:
            formula = db.loc[target_id, 'formula']
            nat = db.loc[target_id, 'num_atoms']
            med, mean, sd_w, mad_w, netsize = np.nan, np.nan, np.nan, np.nan, np.nan
        else:
            df_out = df_out.drop(columns=['d_rxn'])
            mean, med = wmean_wmed(df_out)
            w, v = df_out['weight'], df_out['value']
            sum_w = np.sum(w)
            avg_w = np.sum(w * v) / sum_w
            sd_w = np.sqrt(np.sum(w * (v - avg_w) ** 2) / sum_w)
            mad_w = weighted_median(abs(v - med), w)
            netsize = len(df_out.index)
            formula = db.loc[target_id, 'formula']
            nat = db.loc[target_id, 'num_atoms']
            w *= 100 / np.sum(w)  # normalize weights to 100 (percent)
        self.rn_out = df_out
        vals = target_id, formula, nat, med, mean, sd_w, mad_w, netsize, filter_type, preferred_method, fallback_method
        headers = ['mpid', 'formula', 'num_atoms', 'Prediction(RN)[median]', 'Prediction(RN)[mean]', 'MAD', 'SD',
                   'n_reactions', 'mode', 'preferred_method', 'fallback_method']
        data.append(dict(zip(headers, vals)))
        if verbose:
            print(df_out.to_markdown())
        if save:
            if not os.path.exists(save_path):
                os.makedirs(save_path)
            pd.to_pickle(df_out, os.path.join(save_path, target_id + '.pkl.gz'))
        return data
