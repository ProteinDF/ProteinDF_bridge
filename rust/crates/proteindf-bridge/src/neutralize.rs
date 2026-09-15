// SPDX-FileCopyrightText: The ProteinDF development team
// SPDX-License-Identifier: GPL-3.0-or-later

use crate::atom_group::AtomGroup;
use crate::error::Result;
use crate::ion_pair::IonPair;
use crate::modeling::Modeling;

/// Neutralizes protein charges by adding counter-ions (Na+ / Cl-) to charged residues and termini.
/// Corresponding to `proteindf_bridge.neutralize.Neutralize`.
#[derive(Debug, Clone)]
pub struct Neutralize {
    neutral_obj: AtomGroup,
}

impl Neutralize {
    /// Creates a new `Neutralize` processor and neutralizes the given protein model.
    /// Uses the default embedded reference structures in `Modeling`.
    pub fn new(protein: &AtomGroup) -> Result<Self> {
        let modeling = Modeling::new()?;
        Self::new_with_modeling(protein, &modeling)
    }

    /// Creates a new `Neutralize` processor with a specific `Modeling` instance.
    pub fn new_with_modeling(protein: &AtomGroup, modeling: &Modeling) -> Result<Self> {
        let neutral_obj = Self::neutralize_protein(protein, modeling)?;
        Ok(Self { neutral_obj })
    }

    /// Returns a reference to the neutralized `AtomGroup`.
    pub fn neutralized(&self) -> &AtomGroup {
        &self.neutral_obj
    }

    /// Consumes this processor and returns the neutralized `AtomGroup`.
    pub fn into_neutralized(self) -> AtomGroup {
        self.neutral_obj
    }

    /// Divides an `AtomGroup` path into `(chain_name, res_name)` matching Python's `_divide_path`.
    pub fn divide_path(path: &str) -> (String, String) {
        let parts = AtomGroup::divide_path(path);
        if parts.len() >= 2 {
            (
                parts[parts.len() - 2].clone(),
                parts[parts.len() - 1].clone(),
            )
        } else if parts.len() == 1 {
            (String::new(), parts[0].clone())
        } else {
            (String::new(), String::new())
        }
    }

    /// Computes the exempt list for a model using `IonPair` detection.
    ///
    /// NOTE on divergence / Python parity:
    /// In original Python `neutralize.py`, `_exempt_list()` is defined and tested,
    /// but the call in `_neutralize()` is commented out (`exempt_list = [] # self._exempt_list()`).
    /// Therefore, the exempt list is practically dead code during actual neutralization.
    /// This method is provided to faithfully reproduce Python's behavior and API for testing.
    pub fn exempt_list(model: &AtomGroup) -> Vec<(String, String, String)> {
        let ip = IonPair::new(model);
        let ionpairs = ip.get_ion_pairs();

        let mut exempt = Vec::new();
        for record in ionpairs {
            let (anion_chain, anion_res) = Self::divide_path(&record.anion_path);
            let (cation_chain, cation_res) = Self::divide_path(&record.cation_path);
            exempt.push((anion_chain, anion_res, record.anion_type));
            exempt.push((cation_chain, cation_res, record.cation_type));
        }
        exempt
    }

    /// Neutralizes the protein model non-destructively, returning a new `AtomGroup`.
    fn neutralize_protein(protein: &AtomGroup, modeling: &Modeling) -> Result<AtomGroup> {
        let mut result = protein.clone();
        // In Python neutralize.py: exempt_list = [] # self._exempt_list()
        // The exempt list mechanism is intentionally empty here matching Python.
        let exempt_list: Vec<(String, String, String)> = Vec::new();

        for (model_name, model) in result.groups_mut() {
            log::info!("model: {}", model_name);
            for (chain_name, chain) in model.groups_mut() {
                log::info!("chain: {}", chain_name);
                for (resid, res) in chain.groups_mut() {
                    let resname = res.name.clone();
                    log::info!("res: {}-{}", resid, resname);

                    if res.has_atom("H3") {
                        if !exempt_list
                            .iter()
                            .any(|(c, r, t)| c == chain_name && r == resid && t == "NTM")
                        {
                            let ag = modeling.neutralize_Nterm(res)?;
                            log::info!("add ion for N-term: {:?}", ag);
                            Self::add_ions(res, &ag);
                        } else {
                            log::info!("exempt adding ion: {}/{} Nterm", chain_name, resname);
                        }
                    }

                    if res.has_atom("OXT") {
                        if !exempt_list
                            .iter()
                            .any(|(c, r, t)| c == chain_name && r == resid && t == "CTM")
                        {
                            let ag = modeling.neutralize_Cterm(res)?;
                            log::info!("add ion for C-term: {:?}", ag);
                            Self::add_ions(res, &ag);
                        } else {
                            log::info!("exempt adding ion: {}/{} Cterm", chain_name, resname);
                        }
                    }

                    match resname.as_str() {
                        "GLU" => {
                            if !exempt_list
                                .iter()
                                .any(|(c, r, t)| c == chain_name && r == resid && t == "GLU")
                            {
                                let ag = modeling.neutralize_GLU(res)?;
                                log::info!("add ion for GLU({}): {:?}", resid, ag);
                                Self::add_ions(res, &ag);
                            } else {
                                log::info!("exempt adding ion: {}/{} GLU", chain_name, resname);
                            }
                        }
                        "ASP" => {
                            if !exempt_list
                                .iter()
                                .any(|(c, r, t)| c == chain_name && r == resid && t == "ASP")
                            {
                                let ag = modeling.neutralize_ASP(res)?;
                                log::info!("add ion for ASP({}): {:?}", resid, ag);
                                Self::add_ions(res, &ag);
                            } else {
                                log::info!("exempt adding ion: {}/{} ASP", chain_name, resname);
                            }
                        }
                        "LYS" => {
                            if !exempt_list
                                .iter()
                                .any(|(c, r, t)| c == chain_name && r == resid && t == "LYS")
                            {
                                let ag = modeling.neutralize_LYS(res)?;
                                log::info!("add ion for LYS({}): {:?}", resid, ag);
                                Self::add_ions(res, &ag);
                            } else {
                                log::info!("exempt adding ion: {}/{} LYS", chain_name, resname);
                            }
                        }
                        "ARG" => {
                            let is_exempt = exempt_list.iter().any(|(c, r, t)| {
                                c == chain_name
                                    && r == resid
                                    && (t == "ARG" || t == "ARG1" || t == "ARG2")
                            });
                            if !is_exempt {
                                let ag = modeling.neutralize_ARG(res, 0)?;
                                log::info!("add ion for ARG({}): {:?}", resid, ag);
                                Self::add_ions(res, &ag);
                            } else {
                                log::info!("exempt adding ion: {}/{} ARG", chain_name, resname);
                            }
                        }
                        "FAD" => {
                            let ag = modeling.neutralize_FAD(res)?;
                            log::info!("add ion for FAD({}): {:?}", resid, ag);
                            Self::add_ions(res, &ag);
                        }
                        _ => {}
                    }
                }
            }
        }

        Ok(result)
    }

    /// Adds counter-ions into `atomgroup`, avoiding name collisions matching Python `_add_ions`.
    fn add_ions(atomgroup: &mut AtomGroup, ions: &AtomGroup) {
        let mut count = 0;
        for (atom_name, atom) in ions.atoms() {
            let mut new_name;
            loop {
                new_name = format!("{}{}", atom_name, count);
                if !atomgroup.has_atomkey(&new_name) {
                    break;
                }
                count += 1;
            }
            atomgroup.set_atom(&new_name, atom.clone());
        }
    }
}
