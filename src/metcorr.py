    def analyze(self, event):
        """process event, return True (go to next module) or False (fail, go to next event)"""
        jets      = Collection(event, self.jetBranchName )
        nJet      = event.nJet
        lowPtJets = Collection(event, "CorrT1METJet" ) if self.isV5NanoAOD else []
        muons     = Collection(event, "Muon" ) # to subtract out of the jets for proper type-1 MET corrections
        if not self.isData:
          genJets   = Collection(event, self.genJetBranchName )
        
        # prepare the low pt jets (they don't have a rawFactor)
        for jet in lowPtJets:
          jet.pt        = jet.rawPt
          jet.rawFactor = 0
          jet.mass      = 0
          # the following dummy values should be removed once the values are kept in nanoAOD
          jet.neEmEF    = 0
          jet.chEmEF    = 0

        if not self.isData:
            self.jetSmearer.setSeed(event)
        
        jets_pt_raw = []
        jets_pt_jer = []
        jets_pt_nom = []

        jets_mass_raw = []
        jets_mass_nom = []
        
        jets_corr_JEC = []
        jets_corr_JER = []
        
        jets_pt_jerUp   = []
        jets_pt_jerDown = []
        jets_pt_jesUp   = {}
        jets_pt_jesDown = {}
        
        jets_mass_jerUp   = []
        jets_mass_jerDown = []
        jets_mass_jesUp   = {}
        jets_mass_jesDown = {}

        for jesUncertainty in self.jesUncertainties:
            jets_pt_jesUp[jesUncertainty]   = []
            jets_pt_jesDown[jesUncertainty] = []
            jets_mass_jesUp[jesUncertainty]   = []
            jets_mass_jesDown[jesUncertainty] = []
        
        met     = Object(event, self.metBranchName)
        rawmet  = Object(event, "RawMET")
        defmet  = Object(event, "MET")

        ( t1met_px,       t1met_py       ) = ( met.pt*math.cos(met.phi), met.pt*math.sin(met.phi) )
        ( def_met_px,     def_met_py     ) = ( defmet.pt*math.cos(defmet.phi),   defmet.pt*math.sin(defmet.phi) )
        ( met_px,         met_py         ) = ( rawmet.pt*math.cos(rawmet.phi), rawmet.pt*math.sin(rawmet.phi) )
        ( met_px_nom,     met_py_nom     ) = ( met_px, met_py )
        ( met_px_jer,     met_py_jer     ) = ( met_px, met_py )
        ( met_px_jerUp,   met_py_jerUp   ) = ( met_px, met_py )
        ( met_px_jerDown, met_py_jerDown ) = ( met_px, met_py )
        ( met_px_jesUp,   met_py_jesUp   ) = ( {}, {} )
        ( met_px_jesDown, met_py_jesDown ) = ( {}, {} )

        for jesUncertainty in self.jesUncertainties:
            met_px_jesUp[jesUncertainty]   = met_px
            met_py_jesUp[jesUncertainty]   = met_py
            met_px_jesDown[jesUncertainty] = met_px
            met_py_jesDown[jesUncertainty] = met_py

        # variables needed for re-applying JECs to 2017 v2 MET
        delta_x_T1Jet, delta_y_T1Jet = 0, 0
        delta_x_rawJet, delta_y_rawJet = 0, 0
        
        rho = getattr(event, self.rhoBranchName)
        
        # match reconstructed jets to generator level ones
        # (needed to evaluate JER scale factors and uncertainties)
        def resolution_matching(jet, genjet):
          '''Helper function to match to gen based on pt difference'''
          params = ROOT.PyJetParametersWrapper()
          params.setJetEta(jet.eta)
          params.setJetPt(jet.pt)
          params.setRho(rho)

          resolution = self.jetSmearer.jer.getResolution(params)

          return abs(jet.pt - genjet.pt) < 3*resolution*jet.pt

        if not self.isData:
          pairs = matchObjectCollection(jets, genJets, dRmax=0.2, presel=resolution_matching)
          lowPtPairs = matchObjectCollection(lowPtJets, genJets, dRmax=0.2, presel=resolution_matching)
          pairs.update(lowPtPairs)

        for iJet, jet in enumerate(itertools.chain(jets, lowPtJets)):
            #jet pt and mass corrections
            jet_pt = jet.pt
            jet_mass = jet.mass
            jet_pt_orig = jet_pt
            rawFactor = jet.rawFactor

            #redo JECs if desired
            if hasattr(jet, "rawFactor"):
                jet_rawpt = jet_pt * (1 - jet.rawFactor)
                jet_rawmass = jet_mass * (1 - jet.rawFactor)
            else:
                jet_rawpt = -1.0 * jet_pt #If factor not present factor will be saved as -1
                jet_rawmass = -1.0 * jet_mass #If factor not present factor will be saved as -1

            (jet_pt, jet_mass)    = self.jetReCalibrator.correct(jet,rho)
            (jet_pt_l1, jet_mass_l1) = self.jetReCalibratorL1.correct(jet,rho)
            jet.pt = jet_pt
            jet.mass = jet_mass

            # Get the JEC factors
            jec   = jet_pt/jet_rawpt
            jecL1 = jet_pt_l1/jet_rawpt
            if self.jetReCalibratorProd:
              jecProd = self.jetReCalibratorProd.correct(jet,rho)[0]/jet_rawpt
              jecL1Prod = self.jetReCalibratorProdL1.correct(jet,rho)[0]/jet_rawpt

            if not self.isData:
              genJet = pairs[jet]
            
            # get the jet for type-1 MET
            newjet = ROOT.TLorentzVector()
            if self.isV5NanoAOD:
                newjet.SetPtEtaPhiM(jet_pt_orig*(1-jet.rawFactor)*(1-jet.muonSubtrFactor), jet.eta, jet.phi, jet.mass )
                muon_pt = jet_pt_orig*(1-jet.rawFactor)*jet.muonSubtrFactor
            else:
                newjet.SetPtEtaPhiM(jet_pt_orig*(1-jet.rawFactor), jet.eta, jet.phi, jet.mass )
                muon_pt = 0
                if hasattr(jet, 'muonIdx1'):
                  if jet.muonIdx1>-1:
                      if muons[jet.muonIdx1].isGlobal:
                        newjet = newjet - muons[jet.muonIdx1].p4()
                        muon_pt += muons[jet.muonIdx1].pt
                  if jet.muonIdx2>-1:
                      if muons[jet.muonIdx2].isGlobal:
                        newjet = newjet - muons[jet.muonIdx2].p4()
                        muon_pt += muons[jet.muonIdx2].pt

            # set the jet pt to the muon subtracted raw pt
            jet.pt = newjet.Pt()
            jet.rawFactor = 0
            # get the proper jet pts for type-1 MET. only correct the non-mu fraction of the jet. if the corrected pt>15, use the corrected jet, otherwise use raw
            jet_pt_noMuL1L2L3 = jet.pt*jec    if jet.pt*jec > self.unclEnThreshold else jet.pt
            jet_pt_noMuL1     = jet.pt*jecL1  if jet.pt*jec > self.unclEnThreshold else jet.pt

            # this step is only needed for v2 MET in 2017 when different JECs are applied compared to the nanoAOD production
            if self.jetReCalibratorProd:
              jet_pt_noMuProdL1L2L3   = jet.pt*jecProd    if jet.pt*jecProd > self.unclEnThreshold else jet.pt
              jet_pt_noMuProdL1       = jet.pt*jecL1Prod  if jet.pt*jecProd > self.unclEnThreshold else jet.pt
            else:
              jet_pt_noMuProdL1L2L3 = jet_pt_noMuL1L2L3
              jet_pt_noMuProdL1     = jet_pt_noMuL1

            ## setting jet back to central values
            jet.pt          = jet_pt
            jet.rawFactor   = rawFactor

            # evaluate JER scale factors and uncertainties
            # (cf. https://twiki.cern.ch/twiki/bin/view/CMS/JetResolution and https://twiki.cern.ch/twiki/bin/view/CMSPublic/WorkBookJetEnergyResolution )
            if not self.isData:
              ( jet_pt_jerNomVal, jet_pt_jerUpVal, jet_pt_jerDownVal ) = self.jetSmearer.getSmearValsPt(jet, genJet, rho)
            else:
              # if you want to do something with JER in data, please add it here.
              ( jet_pt_jerNomVal, jet_pt_jerUpVal, jet_pt_jerDownVal ) = (1,1,1)
            
            # these are the important jet pt values
            #jet_pt_nom      = jet_pt if jet_pt > 0 else 0
            jet_pt_nom      = jet_pt * jet_pt_jerNomVal if self.applySmearing else jet_pt
            jet_pt_L1L2L3   = jet_pt_noMuL1L2L3 + muon_pt
            jet_pt_L1       = jet_pt_noMuL1     + muon_pt

            # not nice, but needed for METv2 in 2017
            jet_pt_prodL1L2L3   = jet_pt_noMuProdL1L2L3 + muon_pt
            jet_pt_prodL1       = jet_pt_noMuProdL1     + muon_pt

            if self.metBranchName == 'METFixEE2017':
                # get the delta for removing L1L2L3-L1 corrected jets (corrected with GT from nanoAOD production!!) in the EE region from the default MET branch.
                if jet_pt_prodL1L2L3 > self.unclEnThreshold and 2.65<abs(jet.eta)<3.14 and jet_rawpt < 50:
                    delta_x_T1Jet  += (jet_pt_prodL1L2L3-jet_pt_prodL1) * math.cos(jet.phi) + jet_rawpt * math.cos(jet.phi)
                    delta_y_T1Jet  += (jet_pt_prodL1L2L3-jet_pt_prodL1) * math.sin(jet.phi) + jet_rawpt * math.sin(jet.phi)

                # get the delta for removing raw jets in the EE region from the raw MET
                if jet_pt_prodL1L2L3 > self.unclEnThreshold and 2.65<abs(jet.eta)<3.14 and jet_rawpt < 50:
                    delta_x_rawJet += jet_rawpt * math.cos(jet.phi)
                    delta_y_rawJet += jet_rawpt * math.sin(jet.phi)



            # don't store the low pt jets in the Jet_pt_nom branch
            if iJet < nJet:
                jets_pt_raw     .append(jet_rawpt)
                jets_pt_nom     .append(jet_pt_nom)
                jets_mass_raw   .append(jet_rawmass)
                jets_corr_JEC   .append(jet_pt/jet_rawpt)
                jets_corr_JER   .append(jet_pt_jerNomVal)  # can be used to undo JER

                # no need to do this for low pt jets
                jet_mass_nom           = jet_pt_jerNomVal*jet_mass if self.applySmearing else jet_mass
                if jet_mass_nom < 0.0:
                    jet_mass_nom *= -1.0
                jets_mass_nom    .append(jet_mass_nom)

            if not self.isData:
              jet_pt_jerUp         = jet_pt_jerUpVal  *jet_pt
              jet_pt_jerDown       = jet_pt_jerDownVal*jet_pt

              # evaluate JES uncertainties
              jet_pt_jesUp     = {}
              jet_pt_jesDown   = {}
              jet_pt_jesUpT1   = {}
              jet_pt_jesDownT1 = {}

              jet_mass_jesUp   = {}
              jet_mass_jesDown = {}

              # don't store the low pt jets in the Jet_pt_nom branch
              if iJet < nJet:
                  jets_pt_jerUp    .append(jet_pt_jerUpVal*jet_pt)
                  jets_pt_jerDown  .append(jet_pt_jerDownVal*jet_pt)
                  jets_mass_jerUp  .append(jet_pt_jerUpVal   *jet_mass)
                  jets_mass_jerDown.append(jet_pt_jerDownVal *jet_mass)
              
              for jesUncertainty in self.jesUncertainties:
                  # (cf. https://twiki.cern.ch/twiki/bin/view/CMSPublic/WorkBookJetEnergyCorrections#JetCorUncertainties )
                  self.jesUncertainty[jesUncertainty].setJetPt(jet_pt_nom)
                  self.jesUncertainty[jesUncertainty].setJetEta(jet.eta)
                  delta = self.jesUncertainty[jesUncertainty].getUncertainty(True)
                  jet_pt_jesUp[jesUncertainty]   = jet_pt_nom*(1. + delta)
                  jet_pt_jesDown[jesUncertainty] = jet_pt_nom*(1. - delta)
                  if iJet < nJet:
                    jets_pt_jesUp[jesUncertainty].append(jet_pt_jesUp[jesUncertainty])
                    jets_pt_jesDown[jesUncertainty].append(jet_pt_jesDown[jesUncertainty])
                    jet_mass_jesUp   [jesUncertainty] = jet_mass_nom*(1. + delta)
                    jet_mass_jesDown [jesUncertainty] = jet_mass_nom*(1. - delta)
                    jets_mass_jesUp  [jesUncertainty].append(jet_mass_jesUp[jesUncertainty])
                    jets_mass_jesDown[jesUncertainty].append(jet_mass_jesDown[jesUncertainty])
                  
                  # redo JES variations for T1 MET
                  self.jesUncertainty[jesUncertainty].setJetPt(jet_pt_L1L2L3)
                  self.jesUncertainty[jesUncertainty].setJetEta(jet.eta)
                  delta = self.jesUncertainty[jesUncertainty].getUncertainty(True)
                  jet_pt_jesUpT1[jesUncertainty]   = jet_pt_L1L2L3*(1. + delta)
                  jet_pt_jesDownT1[jesUncertainty] = jet_pt_L1L2L3*(1. - delta)


            # progate JER and JES corrections and uncertainties to MET
            if jet_pt_L1L2L3 > self.unclEnThreshold and (jet.neEmEF+jet.chEmEF) < 0.9:
                if not ( self.metBranchName == 'METFixEE2017' and 2.65<abs(jet.eta)<3.14 and jet.pt*(1-jet.rawFactor)<50 ): # do not re-correct for jets that aren't included in METv2 recipe
                    jet_cosPhi = math.cos(jet.phi)
                    jet_sinPhi = math.sin(jet.phi)
                    met_px_nom     = met_px_nom     - (jet_pt_L1L2L3  - jet_pt_L1)*jet_cosPhi 
                    met_py_nom     = met_py_nom     - (jet_pt_L1L2L3  - jet_pt_L1)*jet_sinPhi 
                    if not self.isData:
                      met_px_jer     = met_px_jer     - (jet_pt_L1L2L3*jet_pt_jerNomVal   - jet_pt_L1)*jet_cosPhi 
                      met_py_jer     = met_py_jer     - (jet_pt_L1L2L3*jet_pt_jerNomVal   - jet_pt_L1)*jet_sinPhi 
                      met_px_jerUp   = met_px_jerUp   - (jet_pt_L1L2L3*jet_pt_jerUpVal    - jet_pt_L1)*jet_cosPhi 
                      met_py_jerUp   = met_py_jerUp   - (jet_pt_L1L2L3*jet_pt_jerUpVal    - jet_pt_L1)*jet_sinPhi 
                      met_px_jerDown = met_px_jerDown - (jet_pt_L1L2L3*jet_pt_jerDownVal  - jet_pt_L1)*jet_cosPhi 
                      met_py_jerDown = met_py_jerDown - (jet_pt_L1L2L3*jet_pt_jerDownVal  - jet_pt_L1)*jet_sinPhi 
                      for jesUncertainty in self.jesUncertainties:
                          met_px_jesUp[jesUncertainty]   = met_px_jesUp[jesUncertainty]   - (jet_pt_jesUpT1[jesUncertainty]   - jet_pt_L1)*jet_cosPhi
                          met_py_jesUp[jesUncertainty]   = met_py_jesUp[jesUncertainty]   - (jet_pt_jesUpT1[jesUncertainty]   - jet_pt_L1)*jet_sinPhi
                          met_px_jesDown[jesUncertainty] = met_px_jesDown[jesUncertainty] - (jet_pt_jesDownT1[jesUncertainty] - jet_pt_L1)*jet_cosPhi
                          met_py_jesDown[jesUncertainty] = met_py_jesDown[jesUncertainty] - (jet_pt_jesDownT1[jesUncertainty] - jet_pt_L1)*jet_sinPhi


        # propagate "unclustered energy" uncertainty to MET
        if self.metBranchName == 'METFixEE2017':
            # Remove the L1L2L3-L1 corrected jets in the EE region from the default MET branch
            def_met_px += delta_x_T1Jet
            def_met_py += delta_y_T1Jet

            # get unclustered energy part that is removed in the v2 recipe
            met_unclEE_x = def_met_px - t1met_px
            met_unclEE_y = def_met_py - t1met_py

            # finalize the v2 recipe for the rawMET by removing the unclustered part in the EE region
            met_px_nom += delta_x_rawJet - met_unclEE_x 
            met_py_nom += delta_y_rawJet - met_unclEE_y
            
            if not self.isData:
              # apply v2 recipe correction factor also to JER central value and variations
              met_px_jer += delta_x_rawJet - met_unclEE_x
              met_py_jer += delta_y_rawJet - met_unclEE_y
              met_px_jerUp += delta_x_rawJet - met_unclEE_x
              met_py_jerUp += delta_y_rawJet - met_unclEE_y
              met_px_jerDown += delta_x_rawJet - met_unclEE_x
              met_py_jerDown += delta_y_rawJet - met_unclEE_y
              for jesUncertainty in self.jesUncertainties:
                  met_px_jesUp[jesUncertainty] += delta_x_rawJet - met_unclEE_x
                  met_py_jesUp[jesUncertainty] += delta_y_rawJet - met_unclEE_y
                  met_px_jesDown[jesUncertainty] += delta_x_rawJet - met_unclEE_x
                  met_py_jesDown[jesUncertainty] += delta_y_rawJet - met_unclEE_y


        if not self.isData:
          ( met_px_unclEnUp,   met_py_unclEnUp   ) = ( met_px_nom, met_py_nom )
          ( met_px_unclEnDown, met_py_unclEnDown ) = ( met_px_nom, met_py_nom )
          met_deltaPx_unclEn = getattr(event, self.metBranchName + "_MetUnclustEnUpDeltaX")
          met_deltaPy_unclEn = getattr(event, self.metBranchName + "_MetUnclustEnUpDeltaY")
          met_px_unclEnUp    = met_px_unclEnUp   + met_deltaPx_unclEn
          met_py_unclEnUp    = met_py_unclEnUp   + met_deltaPy_unclEn
          met_px_unclEnDown  = met_px_unclEnDown - met_deltaPx_unclEn
          met_py_unclEnDown  = met_py_unclEnDown - met_deltaPy_unclEn