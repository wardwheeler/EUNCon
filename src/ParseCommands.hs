{- |
Module      :  ParseCommands.hs
Description :  Progam to calcualte Edge Union Network ala Miyagi and Wheeler 2019
               input graphviz dot files and newick
Copyright   :  (c) 2020 Ward C. Wheeler, Division of Invertebrate Zoology, AMNH. All rights reserved.
License     :

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice, this
   list of conditions and the following disclaimer.
2. Redistributions in binary form must reproduce the above copyright notice,
   this list of conditions and the following disclaimer in the documentation
   and/or other materials provided with the distribution.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR
ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
(INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

The views and conclusions contained in the software and documentation are those
of the authors and should not be interpreted as representing official policies,
either expressed or implied, of the FreeBSD Project.

Maintainer  :  Ward Wheeler <wheeler@amnh.org>
Stability   :  unstable
Portability :  portable (I hope)

Todo:
  from input file?

-}

module ParseCommands (processCommands) where

import qualified Data.Text.Lazy         as T
import           Debug.Trace
import Debug.Trace

-- | processCommands takes a list of strings and returns values of commands for proram execution
-- including defaults
processCommands :: [String] -> (String, String, Int, String, String, [String])
processCommands inList =
    if null inList then error ("No input parameters.\nParameters that can be set:"
        ++ "\n\tMethod=(eun|strict|majority|Adams) "
        ++ "\n\tThreshold=0-100 "
        ++ "\n\tOutFormat=Dot|FENewick"
        ++ "\n\tOutFile=filename"
        ++ "\n\tInput files (may including wildcards) without preceeding \"option=\""
        ++ "\n\tNeed at least a single input graph file (and at least two input graphs)."
        ++ "\n\tDefault values reconcile=EUN, compare=combinable threeshold=0, outformat=dot, outfile=euncon.out\n\n")
    else
        let inTextList = fmap T.pack inList
            inTextListLC = fmap T.toLower inTextList
            inputFileList = getInputFileNames inTextList
            method = getMethod inTextListLC
            compareMethod = getCompareMethod inTextListLC
            threshold = getThreshold inTextListLC
            outFormat = getOutputFormat inTextListLC
            outFile =  getOutputFileName (zip inTextListLC inTextList)
        in
        trace ("\nInput arguments: " ++ (show inList) ++ "\nProgram options: " ++ show (method, compareMethod, threshold, outFormat, outFile, inputFileList))
        (method, compareMethod, threshold, outFormat, outFile, inputFileList)

-- | getInputFileNames returns names not including a parameter '='
getInputFileNames :: [T.Text] -> [String]
getInputFileNames inTextList = T.unpack <$> filter (T.all (/= '=')) inTextList

-- | getMethod returns method value or dedfault otherwise
-- assumes in lower case
getMethod :: [T.Text] -> String
getMethod inTextList =
    if null inTextList then trace ("Warning: No reconcile specified defaulting to \'eun\'") "eun"
    else
        let firstCommand = T.takeWhile (/= '=') $ head inTextList
            firstOption = T.tail $ T.dropWhile (/= '=') $ head inTextList
        in
        if firstCommand == T.pack "reconcile" then 
            let option = T.unpack firstOption
            in 
            if option == "eun" then "eun"
            else if option == "majority" then "majority"
            else if option == "strict" then "strict"
            else error ("Reconcile option \'" ++ option ++ "\' not recognized")
        else getMethod (tail inTextList)

-- | getCompareMethod returns compareMethod value or dedfault otherwise
-- assumes in lower case
getCompareMethod :: [T.Text] -> String
getCompareMethod inTextList =

    if null inTextList then trace ("Warning: No compare specified defaulting to \'combinable\'") "combinable"
    else
        let firstCommand = T.takeWhile (/= '=') $ head inTextList
            firstOption = T.tail $ T.dropWhile (/= '=') $ head inTextList
        in
        if firstCommand == T.pack "compare" then 
            let option = T.unpack firstOption
            in 
            if option == "combinable" then "combinable"
            else if option == "identity" then"identity"
            else error ("Compare option \'" ++ option ++ "\' not recognized")
        else getCompareMethod (tail inTextList)

-- | getThreshold returns threshold value or default otherwise
-- assumes in lower case
getThreshold :: [T.Text] -> Int
getThreshold inTextList =
    if null inTextList then trace ("Warning: No threshold specified defaulting to \'0\'")  0 :: Int
    else
        let firstCommand = T.takeWhile (/= '=') $ head inTextList
            firstOption = T.tail $ T.dropWhile (/= '=') $ head inTextList
        in
        if firstCommand == T.pack "threshold" then read (T.unpack firstOption) :: Int
        else getThreshold (tail inTextList)

-- | getOutputFormat returns output file format or default otherwise
-- assumes in lower case
getOutputFormat :: [T.Text] -> String
getOutputFormat inTextList =
    if null inTextList then trace ("Warning: No output format specified defaulting to \'dot\'") "dot"
    else
        let firstCommand = T.takeWhile (/= '=') $ head inTextList
            firstOption = T.tail $ T.dropWhile (/= '=') $ head inTextList
        in
        if firstCommand == T.pack "outformat" then T.unpack firstOption
        else getOutputFormat (tail inTextList)

-- | getOutputFileName returns output file name or default otherwise
-- assumes in lower case for command, uses pair so no case convewrsino in files name
getOutputFileName :: [(T.Text, T.Text)] -> String
getOutputFileName inTextPairList =
    if null inTextPairList then trace ("Warning: No output file name specified defaulting to \'euncon.out\'") "euncon.out"
    else
        let (textListLC, textList) = head inTextPairList
            firstCommand = T.takeWhile (/= '=') textListLC
            firstOption = T.tail $ T.dropWhile (/= '=') textList
        in
        if firstCommand == T.pack "outfile" then T.unpack firstOption
        else getOutputFileName (tail inTextPairList)

