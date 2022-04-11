#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri May 14 09:34:26 2021

@author: mattia
"""


class EmptyRanking:
    def __init__(self):
        self.source : str
        self.ranking = {
            "H200": [0, "Unstable explosive"],
            "H201": [0, "Explosive; mass explosion hazard"],
            "H202": [0, "Explosive; severe projection hazard"],
            "H203": [0, "Explosive; fire, blast or projection hazard"],
            "H204": [0, "Fire or projection hazard"],
            "H205": [0, "May mass explode in fire"],
            "H206": [
                0,
                "Fire, blast or projection hazard; increased risk of explosion if desensitizing agent is reduced",
            ],
            "H207": [
                0,
                "Fire or projection hazard; increased risk of explosion if desensitizing agent is reduced",
            ],
            "H208": [
                0,
                "Fire hazard; increased risk of explosion if desensitizing agent is reduced",
            ],
            "H220": [0, "Extremely flammable gas"],
            "H221": [0, "Flammable gas"],
            "H222": [0, "Extremely flammable aerosol"],
            "H223": [0, "Flammable aerosol"],
            "H224": [0, "Extremely flammable liquid and vapour"],
            "H225": [0, "Highly flammable liquid and vapour"],
            "H226": [0, "Flammable liquid and vapour"],
            "H227": [0, "Combustible liquid"],
            "H228": [0, "Flammable solid"],
            "H229": [0, "Pressurized container: may burst if heated"],
            "H230": [0, "May react explosively even in the absence of air"],
            "H231": [
                0,
                "May react explosively even in the absence of air at elevated pressure and/or temperature",
            ],
            "H232": [0, "May ignite spontaneously if exposed to air"],
            "H240": [0, "Heating may cause an explosion"],
            "H241": [0, "Heating may cause a fire or explosion"],
            "H242": [0, "Heating may cause a fire"],
            "H250": [0, "Catches fire spontaneously if exposed to air"],
            "H251": [0, "Self-heating; may catch fire"],
            "H252": [0, "Self-heating in large quantities; may catch fire"],
            "H260": [
                0,
                "In contact with water releases flammable gases which may ignite spontaneously",
            ],
            "H261": [0, "In contact with water releases flammable gas"],
            "H270": [0, "May cause or intensify fire; oxidizer"],
            "H271": [0, "May cause fire or explosion; strong oxidizer"],
            "H272": [0, "May intensify fire; oxidizer"],
            "H280": [0, "Contains gas under pressure; may explode if heated"],
            "H281": [
                0,
                "Contains refrigerated gas; may cause cryogenic burns or injury",
            ],
            "H282": [
                0,
                "Extremely flammable chemical under pressure: may explode if heated",
            ],
            "H283": [0, "Flammable chemical under pressure: may explode if heated"],
            "H284": [0, "Chemical under pressure: may explode if heated"],
            "H290": [0, "May be corrosive to metals"],
            "H300": [0, "Fatal if swallowed"],
            "H301": [0, "Toxic if swallowed"],
            "H302": [0, "Harmful if swallowed"],
            "H303": [0, "May be harmful if swallowed"],
            "H304": [0, "May be fatal if swallowed and enters airways"],
            "H305": [0, "May be harmful if swallowed and enters airways"],
            "H310": [0, "Fatal in contact with skin"],
            "H311": [0, "Toxic in contact with skin"],
            "H312": [0, "Harmful in contact with skin"],
            "H313": [0, "May be harmful in contact with skin"],
            "H314": [0, "Causes severe skin burns and eye damage"],
            "H315": [0, "Causes skin irritation"],
            "H316": [0, "Causes mild skin irritation"],
            "H317": [0, "May cause an allergic skin reaction"],
            "H318": [0, "Causes serious eye damage"],
            "H319": [0, "Causes serious eye irritation"],
            "H320": [0, "Causes eye irritation"],
            "H330": [0, "Fatal if inhaled"],
            "H331": [0, "Toxic if inhaled"],
            "H332": [0, "Harmful if inhaled"],
            "H333": [0, "May be harmful if inhaled"],
            "H334": [
                0,
                "May cause allergy or asthma symptoms or breathing difficulties if inhaled",
            ],
            "H335": [0, "May cause respiratory irritation"],
            "H336": [0, "May cause drowsiness or dizziness"],
            "H340": [0, "May cause genetic defects"],
            "H341": [0, "Suspected of causing genetic defects"],
            "H350": [0, "May cause cancer"],
            "H350i": [0, "May cause cancer by inhalation"],
            "H351": [0, "Suspected of causing cancer"],
            "H360": [0, "May damage fertility or the unborn child"],
            "H361": [0, "Suspected of damaging fertility or the unborn child"],
            "H361d": [0, "Suspected of damaging the unborn child"],
            "H361D": [0, "May damage the unborn child"],
            "H361f": [0, "Suspected of damaging fertility"],
            "H361F": [0, "May damage fertility"],
            "H362": [0, "May cause harm to breast-fed children"],
            "H370": [0, "Causes damage to organs"],
            "H371": [0, "May cause damage to organs"],
            "H372": [
                0,
                "Causes damage to organs through prolonged or repeated exposure",
            ],
            "H373": [
                0,
                "May cause damage to organs through prolonged or repeated exposure",
            ],
            "H400": [0, "Very toxic to aquatic life"],
            "H401": [0, "Toxic to aquatic life"],
            "H402": [0, "Harmful to aquatic life"],
            "H410": [0, "Very toxic to aquatic life with long-lasting effects"],
            "H411": [0, "Toxic to aquatic life with long-lasting effects"],
            "H412": [0, "Harmful to aquatic life with long-lasting effects"],
            "H413": [0, "May cause long-lasting harmful effects to aquatic life"],
            "H420": [
                0,
                "Harms public health and the environment by destroying ozone in the upper atmosphere",
            ],
        }

    def evaluate(self, references):
        compound_keys = list(references.keys())

        scores = {}

        for compound in compound_keys:
            scores[compound] = self._rank(references[compound])

        return scores

    def _rank(self, reference):
        ref_keys = list(reference.keys())

        # if no hazard phrase is present for the compound, return None
        # -1 indicates the compound was not found on PUBCHEM
        # -2 indicates the compound was found but not GHS info is present
        if (ref_keys == [-1]):
            return -1
        if (ref_keys == [-2]):
            return -2

        # If the ref is empty due to cleaning non-ECHA sources
        if ref_keys == []:
            return None

        # Return sum of rankings of the hazard phrases if an ECHA source if present
        for key in ref_keys:
            if reference[key].source_name == self.source:
                hazard_codes = reference[key].hazard_codes
                break
            else:
                # Don't know if this is useful...
                hazard_codes = []
                continue

        # In case no hazard phrase if found, return None
        # P.s.: there should not be any case where hazard_codes is None, but lets not mess up
        if (hazard_codes is None) or (hazard_codes == []):
            return None

        score = 0

        for code in hazard_codes:
            if code is None: # for non hazardous compounds
                score += 0
            else:
                score += self.ranking[code][0]

        return score


class GES001(EmptyRanking):
    def __init__(self):
        super().__init__()
        self.ranking["H200"][0] = 1000
        self.ranking["H201"][0] = 1000
        self.ranking["H202"][0] = 1000
        self.ranking["H203"][0] = 1000
        self.ranking["H204"][0] = 1000
        self.ranking["H205"][0] = 1000
        self.ranking["H206"][0] = 1000
        self.ranking["H207"][0] = 1000
        self.ranking["H208"][0] = 5
        self.ranking["H220"][0] = 1000
        self.ranking["H221"][0] = 1000
        self.ranking["H222"][0] = 1000
        self.ranking["H223"][0] = 5
        self.ranking["H224"][0] = 1000
        self.ranking["H225"][0] = 1000
        self.ranking["H226"][0] = 1000
        self.ranking["H227"][0] = 5
        self.ranking["H228"][0] = 1000
        self.ranking["H229"][0] = 1000
        self.ranking["H230"][0] = 5
        self.ranking["H231"][0] = 1000
        self.ranking["H232"][0] = 1000
        self.ranking["H240"][0] = 1000
        self.ranking["H241"][0] = 1000
        self.ranking["H242"][0] = 1000
        self.ranking["H250"][0] = 1000
        self.ranking["H251"][0] = 1000
        self.ranking["H252"][0] = 1000
        self.ranking["H260"][0] = 1000
        self.ranking["H261"][0] = 1000
        self.ranking["H270"][0] = 1000
        self.ranking["H271"][0] = 1000
        self.ranking["H272"][0] = 5
        self.ranking["H280"][0] = 5
        self.ranking["H281"][0] = 5
        self.ranking["H282"][0] = 1000
        self.ranking["H283"][0] = 1000
        self.ranking["H284"][0] = 5
        self.ranking["H290"][0] = 5
        self.ranking["H300"][0] = 1000
        self.ranking["H301"][0] = 1000
        self.ranking["H302"][0] = 4
        self.ranking["H303"][0] = 3
        self.ranking["H304"][0] = 1000
        self.ranking["H305"][0] = 1000
        self.ranking["H310"][0] = 1000
        self.ranking["H311"][0] = 1000
        self.ranking["H312"][0] = 5
        self.ranking["H313"][0] = 4
        self.ranking["H314"][0] = 1000
        self.ranking["H315"][0] = 5
        self.ranking["H316"][0] = 4
        self.ranking["H317"][0] = 4
        self.ranking["H318"][0] = 1000
        self.ranking["H319"][0] = 5
        self.ranking["H320"][0] = 4
        self.ranking["H330"][0] = 1000
        self.ranking["H331"][0] = 1000
        self.ranking["H332"][0] = 5
        self.ranking["H333"][0] = 4
        self.ranking["H334"][0] = 1000
        self.ranking["H335"][0] = 4
        self.ranking["H336"][0] = 3
        self.ranking["H340"][0] = 1000
        self.ranking["H341"][0] = 1000
        self.ranking["H350"][0] = 1000
        self.ranking["H350i"][0] = 1000
        self.ranking["H351"][0] = 1000
        self.ranking["H360"][0] = 1000
        self.ranking["H361"][0] = 1000
        self.ranking["H361d"][0] = 1000
        self.ranking["H361D"][0] = 1000
        self.ranking["H361f"][0] = 1000
        self.ranking["H361F"][0] = 1000
        self.ranking["H362"][0] = 5
        self.ranking["H370"][0] = 1000
        self.ranking["H371"][0] = 5
        self.ranking["H372"][0] = 5
        self.ranking["H373"][0] = 5
        self.ranking["H400"][0] = 4
        self.ranking["H401"][0] = 4
        self.ranking["H402"][0] = 3
        self.ranking["H410"][0] = 4
        self.ranking["H411"][0] = 4
        self.ranking["H412"][0] = 3
        self.ranking["H413"][0] = 3
        self.ranking["H420"][0] = 4
        
class GES002(EmptyRanking):
    def __init__(self):
        super().__init__()
        self.source = 'European Chemicals Agency (ECHA)'
        self.ranking["H200"][0] = 1000
        self.ranking["H201"][0] = 1000
        self.ranking["H202"][0] = 1000
        self.ranking["H203"][0] = 1000
        self.ranking["H204"][0] = 1000
        self.ranking["H205"][0] = 1000
        self.ranking["H206"][0] = 1000
        self.ranking["H207"][0] = 1000
        self.ranking["H208"][0] = 5
        self.ranking["H220"][0] = 5 # e.g. hydrogen
        self.ranking["H221"][0] = 5
        self.ranking["H222"][0] = 1000
        self.ranking["H223"][0] = 5
        self.ranking["H224"][0] = 1000
        self.ranking["H225"][0] = 1000
        self.ranking["H226"][0] = 1000
        self.ranking["H227"][0] = 5
        self.ranking["H228"][0] = 1000
        self.ranking["H229"][0] = 3 
        self.ranking["H230"][0] = 5
        self.ranking["H231"][0] = 5
        self.ranking["H232"][0] = 1000
        self.ranking["H240"][0] = 1000
        self.ranking["H241"][0] = 1000
        self.ranking["H242"][0] = 5
        self.ranking["H250"][0] = 1000
        self.ranking["H251"][0] = 1000
        self.ranking["H252"][0] = 1000
        self.ranking["H260"][0] = 1000
        self.ranking["H261"][0] = 1000
        self.ranking["H270"][0] = 5
        self.ranking["H271"][0] = 1000
        self.ranking["H272"][0] = 5
        self.ranking["H280"][0] = 5
        self.ranking["H281"][0] = 5
        self.ranking["H282"][0] = 5
        self.ranking["H283"][0] = 5
        self.ranking["H284"][0] = 5
        self.ranking["H290"][0] = 3
        self.ranking["H300"][0] = 1000
        self.ranking["H301"][0] = 1000
        self.ranking["H302"][0] = 4
        self.ranking["H303"][0] = 3
        self.ranking["H304"][0] = 1000
        self.ranking["H305"][0] = 1000
        self.ranking["H310"][0] = 1000
        self.ranking["H311"][0] = 1000
        self.ranking["H312"][0] = 5
        self.ranking["H313"][0] = 4
        self.ranking["H314"][0] = 5 # like NaOH or H2SO4
        self.ranking["H315"][0] = 5
        self.ranking["H316"][0] = 4
        self.ranking["H317"][0] = 4
        self.ranking["H318"][0] = 5
        self.ranking["H319"][0] = 4
        self.ranking["H320"][0] = 3
        self.ranking["H330"][0] = 1000
        self.ranking["H331"][0] = 1000
        self.ranking["H332"][0] = 5
        self.ranking["H333"][0] = 4
        self.ranking["H334"][0] = 5
        self.ranking["H335"][0] = 4
        self.ranking["H336"][0] = 3
        self.ranking["H340"][0] = 1000
        self.ranking["H341"][0] = 1000
        self.ranking["H350"][0] = 1000
        self.ranking["H350i"][0] = 1000
        self.ranking["H351"][0] = 1000
        self.ranking["H360"][0] = 1000
        self.ranking["H361"][0] = 1000
        self.ranking["H361d"][0] = 1000
        self.ranking["H361D"][0] = 1000
        self.ranking["H361f"][0] = 1000
        self.ranking["H361F"][0] = 1000
        self.ranking["H362"][0] = 1000
        self.ranking["H370"][0] = 1000
        self.ranking["H371"][0] = 5
        self.ranking["H372"][0] = 5
        self.ranking["H373"][0] = 5
        self.ranking["H400"][0] = 3
        self.ranking["H401"][0] = 3
        self.ranking["H402"][0] = 2
        self.ranking["H410"][0] = 3
        self.ranking["H411"][0] = 3
        self.ranking["H412"][0] = 2
        self.ranking["H413"][0] = 2
        self.ranking["H420"][0] = 3

GES002_instance = GES002()