#  This file is part of RxnNet.
#  Copyright © 2024 Technical University of Denmark

import os
import pandas as pd
import numpy as np
from chempy import balance_stoichiometry
from pymatgen.analysis import reaction_calculator
from pymatgen.core import Composition


class GenerateRN:
    """
    Class for generating candidate reactions based on the reaction network (RN) approach.
    """

    def __init__(self, target_id, db=r'..\rxnnet\data\database-20.03.24\datasets\main.pkl.gz', similarity_distance=True, chemistry_balance=True) -> None:
        """

        :param target_id: target compound identifier (str)
        :param db: filename of database containing at minimum the elemental compositions of the compounds (str)
        :param similarity_distance: whether to calculate reaction similarity distance parameter (bool)
        :param chemistry_balance: whether to determine if reaction chemistry is balanced (bool)
        """
        self.target_id = target_id
        self.db = pd.read_pickle(db)
        self.similarity_distance = similarity_distance
        self.chemistry_balance = chemistry_balance

    def find_product_candidates(self):
        """
        Finds candidate products based on least common (relative to database) element.

        :return: list of candidate products (list)
        """
        target_id = self.target_id
        db = self.db

        db = db.to_dict('index')
        stoic = Composition(db[target_id]["Composition"]).elements
        candidates = {}
        for el in stoic:
            candidates[el] = [candidate_id for candidate_id, candidate in db.items() if
                              (el in Composition(candidate['Composition']).elements and candidate_id != target_id)]
        candidates_len = dict((key, len(val)) for key, val in candidates.items())
        candidates_len_sortedkeys = sorted(candidates_len, key=candidates_len.get)
        return candidates[candidates_len_sortedkeys[0]]

    def find_unbalanced_reactions(self):
        """
        Finds unbalanced reactions which satisfy element identity constraints, i.e. same elemental composition on both
        sides of reaction.

        :return: unbalanced reactions (dict)
        """

        target_id = self.target_id
        db = self.db

        # generate dict of elemental compositions (no stoichiometry)
        el_comp = {}
        for mp_id, v in db.to_dict('index').items():
            el_comp[mp_id] = list(Composition(v['Composition']).elements)
        # add blank entry:
        el_comp[''] = []
        el_comp_target = el_comp[target_id]
        # generate initial list of candidate products
        prod_cands_1 = self.find_product_candidates()
        # search for reactions
        reactions, reaction_no = {}, 0
        for prod_id_1 in prod_cands_1:
            el_comp_prod_cand_1 = el_comp[prod_id_1]
            prod_cands_2 = [trial_id for trial_id, trial_comp in el_comp.items() if
                            (set(el_comp_target).issubset(set(el_comp_prod_cand_1 + trial_comp))
                             and trial_id not in (prod_id_1, target_id)) or
                            (set(trial_comp).issubset(el_comp_prod_cand_1) and trial_id not in (prod_id_1, target_id))
                            or (set(el_comp_prod_cand_1).issubset(el_comp_target)
                                and trial_id not in (prod_id_1, target_id))]
            for prod_id_2 in prod_cands_2:
                el_comp_prod_cand_2 = el_comp[prod_id_2]
                el_comp_prod_cand_1_2 = el_comp_prod_cand_1 + el_comp_prod_cand_2
                reactant_cands_1 = [trial_id for trial_id, trial_comp in el_comp.items() if (
                        set(trial_comp).union(set(el_comp_target)) == set(el_comp_prod_cand_1_2)) and trial_id not in (
                                        prod_id_2, prod_id_1, target_id)]
                for reac_id_1 in reactant_cands_1:
                    reacs_id, prods_id = [target_id, reac_id_1], [prod_id_1, prod_id_2]
                    for x_id in [reacs_id, prods_id]:
                        try:
                            x_id.remove('')
                        except ValueError:
                            pass
                    reactions[reaction_no] = {'reactants': {'ids': reacs_id}, 'products': {'ids': prods_id}}
                    reaction_no += 1
        return reactions

    def rn_similarity_distance(self, rn: pd.DataFrame):
        """
        Find reaction similarity distance parameter for all reactions in reaction network dataframe.

        :param rn: processed RN dataframe (pandas Dataframe)
        :return: list of reaction similarity distance parameters (d_rxn) for every reaction in network (list)
        """
        db = self.db

        rn_local = rn.copy()

        def geometrical_similarity(c_id_rxn):
            """

            :param c_id_rxn:
            :return:
            """
            simil_a = []
            for c_id1 in c_id_rxn['reactants']:
                across_a = []
                for c_id2 in c_id_rxn['products']:
                    v_1 = np.array(db.loc[c_id1, "structural_fingerprint"])
                    v_2 = np.array(db.loc[c_id2, "structural_fingerprint"])
                    simil = np.linalg.norm(v_2 - v_1)
                    across_a.append(simil)
                simil_a.append(across_a)
            return simil_a

        rn_local['similarity_distance'] = [geometrical_similarity(c_id_rxn) for c_id_rxn in
                                           rn_local['reaction_mpids']]
        return [np.mean(np.array(simil_a)) for simil_a in rn_local['similarity_distance']]

    def rn_balanced_chemistry(self, rn: pd.DataFrame):
        """
        Determine whether reactions are balanced with respect to chemistry (via their labels)

        :param rn: processed RN dataframe (pandas Dataframe)
        :return: list of chemistry labels and list of bools representing whether reaction is balanced (tuple)
        """
        db = self.db

        rn_local = rn.copy()
        reaction_partners = rn_local['reaction_mpids']
        mask = []
        reaction_chemistry = []

        def get_chemistry_specie(rxn, rxnside=''):
            """

            :param rxn:
            :param rxnside:
            :return:
            """
            at_list = []
            for specie in rxn[rxnside]:
                if specie is not None:
                    at = db.loc[specie, 'chemistry_type']
                    at_list.extend(at)
            return at_list

        for rxn in reaction_partners:
            reactants, products = [get_chemistry_specie(rxn, rxnside) for rxnside in ['reactants', 'products']]
            reactants.sort()
            products.sort()
            mask_value = set(reactants) == set(products)  # relaxed equivalence
            mask.append(mask_value)
            reaction_chemistry.append({'reactants': reactants, 'products': products})
        return reaction_chemistry, mask

    def balanced_reactions(self, save=False, save_path='user_reactions') -> pd.DataFrame:
        """
        Finds and balances reactions given target id and database of candidate compounds.

        :return: reaction network dataframe containing balanced reactions (pandas Dataframe)
        """

        target_id = self.target_id
        db = self.db
        print('Finding candidate unbalanced reactions...')
        reactions = self.find_unbalanced_reactions()
        rxn_no, data = 0, []
        reactions_clean = []
        print('Balancing reactions, determining whether they are unique and discarding ones that cannot be balanced...')
        for rxn in reactions.keys():
            reacs_in, prods_in = reactions[rxn]['reactants']['ids'], reactions[rxn]['products']['ids']
            reacs_in_comp, prods_in_comp = (
                dict([(c_id, Composition(db.loc[c_id, 'Composition'])) for c_id in reacs_in]),
                dict([(c_id, Composition(db.loc[c_id, 'Composition'])) for c_id in prods_in]))
            # try to balance reaction
            try:
                rxn_bal = reaction_calculator.Reaction(list(reacs_in_comp.values()), list(prods_in_comp.values()))
                rxn_bal_id = sorted(rxn_bal.normalized_repr.split(' '))
                balanced = True
            except reaction_calculator.ReactionError:
                rxn_bal = None
                rxn_bal_id = None
                balanced = False
                pass
            if balanced and rxn_bal_id not in reactions_clean:
                species_comp = reacs_in_comp | prods_in_comp
                reactions_clean.append(rxn_bal_id)
                rxn_no += 1
                stoichiometry = dict(
                    [(reacs_in[idx], v) for idx, v in enumerate(rxn_bal.coeffs[:len(reacs_in_comp)])]) | dict(
                    [(prods_in[idx], v) for idx, v in enumerate(rxn_bal.coeffs[len(reacs_in_comp):])])
                reacs, prods = (dict([(k, v) for k, v in stoichiometry.items() if -1e5 < float(round(v, 5)) < 0]),
                                dict([(k, v) for k, v in stoichiometry.items() if 1e5 > float(round(v, 5)) > 0]))
                if target_id in prods:
                    reacs, prods = dict([(k, -v) for k, v in prods.items()]), dict([(k, -v) for k, v in reacs.items()])
                if target_id in reacs:
                    try:
                        # check for and skip ambiguous/underdetermined reactions
                        balance_stoichiometry({species_comp[k].to_pretty_string() for k, v in reacs.items()},
                                              {species_comp[k].to_pretty_string() for k, v in prods.items()},
                                              underdetermined=False)
                        ambiguous = False
                    except ValueError:
                        ambiguous = True
                        pass
                    if not ambiguous:
                        # reformat
                        rxn_bal_out = reaction_calculator.BalancedReaction(
                            dict((species_comp[k], float(abs(v))) for k, v in reacs.items()),
                            dict((species_comp[k], float(abs(v))) for k, v in prods.items()))
                        rxn_bal_out.normalize_to(species_comp[target_id])
                        stoichiometry = {
                            'reactants': tuple([rxn_bal_out.get_coeff(comp) for comp in rxn_bal_out.reactants]),
                            'products': tuple([rxn_bal_out.get_coeff(comp) for comp in rxn_bal_out.products])}
                        stoichiometry_at = {'reactants': tuple(
                            [rxn_bal_out.get_coeff(comp) * comp.num_atoms for comp in rxn_bal_out.reactants]),
                            'products': tuple([rxn_bal_out.get_coeff(comp) * comp.num_atoms for comp in
                                               rxn_bal_out.products])}
                        ids = {'reactants': tuple(
                            [list(species_comp.keys())[list(species_comp.values()).index(comp)] for comp in
                             rxn_bal_out.reactants]), 'products': tuple(
                            [list(species_comp.keys())[list(species_comp.values()).index(comp)] for comp in
                             rxn_bal_out.products])}
                        headers = ('reaction', 'stoichiometry (atoms)', 'reaction_mpids', 'BalancedReaction')
                        vals = (rxn_bal_out, stoichiometry_at, ids, rxn_bal_out.as_dict())
                        data.append(dict(zip(headers, vals)))
        rn = pd.DataFrame(data)
        if self.similarity_distance:
            rn['d_rxn'] = self.rn_similarity_distance(rn)
        if self.chemistry_balance:
            rn['reaction_chemistry'], rn['balanced_chemistry'] = self.rn_balanced_chemistry(rn)
        rn = rn.set_index(np.arange(1, len(rn.index) + 1, 1))
        rn.index.names = ['reaction number']
        if save:
            if not os.path.exists(save_path):
                os.makedirs(save_path)
            pd.to_pickle(rn, os.path.join(save_path, target_id + '.pkl.gz'))
        return rn
