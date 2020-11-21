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
  add edit distance to commands for error

-}

module ParseCommands (processCommands) where

import qualified Data.Text.Lazy         as T
import           Debug.Trace
import           Data.Array

-- | function for first element of triple
fst3 :: (a,b,c) -> a
fst3 (d,_,_) = d

-- | editDistance is a naive edit distance between two lists
-- takes two  lists and returns edit distance
--- from  https://wiki.haskell.org/Edit_distance
editDistance :: Eq a => [a] -> [a] -> Int
editDistance xs ys = table ! (m,n)
    where
    (m,n) = (length xs, length ys)
    x     = array (1,m) (zip [1..] xs)
    y     = array (1,n) (zip [1..] ys)

    table :: Array (Int,Int) Int
    table = array bnds [(ij, dist ij) | ij <- range bnds]
    bnds  = ((0,0),(m,n))

    dist (0,j) = j
    dist (i,0) = i
    dist (i,j) = minimum [table ! (i-1,j) + 1, table ! (i,j-1) + 1,
        if x ! i == y ! j then table ! (i-1,j-1) else 1 + table ! (i-1,j-1)]

-- | allowedCommandList list of allowable commands
allowedCommandList :: [String]
allowedCommandList = ["reconcile", "compare", "threshold", "outformat", "outfile", "connect"]

-- | getBestMatch compares input to allowable commands and checks if in list and if not outputs
-- closest match
getBestMatch :: (Int, String) -> [String] -> String -> (Int, String)
getBestMatch curBest@(minDist, _) allowedStrings inString =
    if null allowedStrings then curBest
    else
        let candidate =  head allowedStrings
            candidateEditCost = editDistance candidate inString
        in
        if candidateEditCost == 0 then (0, candidate)
        else if candidateEditCost < minDist then getBestMatch (candidateEditCost, candidate) (tail allowedStrings) inString
        else getBestMatch curBest (tail allowedStrings) inString

-- | getCommandErrorString takes list of non zero edits to allowed commands and reurns meaningful error string
getCommandErrorString :: [(Int, String, String)] -> String
getCommandErrorString noMatchList =
    if null noMatchList then ""
    else
        let (_, firstCommand, firstMatch) = head noMatchList
            firstError = "\tBy \'" ++ firstCommand ++ "\' did you mean \'" ++ firstMatch ++ "\'?\n"
        in
        firstError ++ getCommandErrorString (tail noMatchList)

-- | processCommands takes a list of strings and returns values of commands for proram execution
-- including defaults
-- checks commands for misspellings
processCommands :: [String] -> (String, String, Int, Bool, String, String, [String])
processCommands inList =
    if null inList then error ("\n\nError--No input parameters.\nParameters that can be set:"
        ++ "\n\tReconcile=(eun|cun|strict|majority|Adams) "
        ++ "\n\tCompare=combinable|identity "
        ++ "\n\tThreshold=0-100 "
        ++ "\n\tOutFormat=Dot|FENewick"
        ++ "\n\tConnect=True|False"
        ++ "\n\tOutFile=filename"
        ++ "\n\tInput files (may include wildcards) without preceeding \"option=\""
        ++ "\n\tRequires at least a single input graph file (and at least two input graphs)."
        ++ "\n\tDefault values reconcile=EUN, compare=combinable threeshold=0, outformat=dot, outfile=euncon.out\n\n")
    else
        let inTextList = fmap T.pack inList
            inTextListLC = fmap T.toLower inTextList
            commandList = filter (T.any (== '=')) inTextListLC
            stringCommands = fmap (T.unpack . T.takeWhile (/= '=')) commandList
            (editCostList, matchList) = unzip $ fmap (getBestMatch (maxBound :: Int ,"no suggestion") allowedCommandList) stringCommands
            commandMatch = zip3 editCostList stringCommands matchList
            notMatchedList = filter ((>0).fst3) commandMatch
            inputFileList = getInputFileNames inTextList
            method = getMethod inTextListLC
            compareMethod = getCompareMethod inTextListLC
            connect = getConnect inTextListLC
            threshold = if method == "cun" then 0 
                        else if method == "strict" then 100 
                        else getThreshold inTextListLC
            outFormat = getOutputFormat inTextListLC
            outFile =  getOutputFileName (zip inTextListLC inTextList)
        in
        if null notMatchedList then
            trace ("\nInput arguments: " ++ show inList ++ "\nProgram options: " ++ show (method, compareMethod, threshold, connect, outFormat, outFile, inputFileList))
            (method, compareMethod, threshold, connect, outFormat, outFile, inputFileList)
        else error ("\n\nError(s) in command specification (case insensitive):\n" ++ getCommandErrorString notMatchedList)


-- | getInputFileNames returns names not including a parameter '='
getInputFileNames :: [T.Text] -> [String]
getInputFileNames inTextList = T.unpack <$> filter (T.all (/= '=')) inTextList

-- | getMethod returns method value or dedfault otherwise
-- assumes in lower case
getMethod :: [T.Text] -> String
getMethod inTextList =
    if null inTextList then trace "Warning: No reconcile specified defaulting to \'eun\'" "eun"
    else
        let firstCommand = T.takeWhile (/= '=') $ head inTextList
            firstOption = T.tail $ T.dropWhile (/= '=') $ head inTextList
        in
        if firstCommand == T.pack "reconcile" then
            let option = T.unpack firstOption
            in
            if option == "eun" then "eun"
            else if option == "cun" then "cun"
            else if option == "majority" then "majority"
            else if option == "strict" then "strict"
            else error ("Reconcile option \'" ++ option ++ "\' not recognized (eun|cun|majority|strict)")
        else getMethod (tail inTextList)

-- | getCompareMethod returns compareMethod value or default otherwise
-- assumes in lower case
getCompareMethod :: [T.Text] -> String
getCompareMethod inTextList =
    if null inTextList then trace "Warning: No compare specified defaulting to \'combinable\'" "combinable"
    else
        let firstCommand = T.takeWhile (/= '=') $ head inTextList
            firstOption = T.tail $ T.dropWhile (/= '=') $ head inTextList
        in
        if firstCommand == T.pack "compare" then
            let option = T.unpack firstOption
            in
            if option == "combinable" then "combinable"
            else if option == "identity" then "identity"
            else error ("Compare option \'" ++ option ++ "\' not recognized (combinable|identity)")
        else getCompareMethod (tail inTextList)

-- | getConect returns connect value or default otherwise (True|False)
-- assumes in lower case
getConnect :: [T.Text] -> Bool
getConnect inTextList =
    if null inTextList then trace "Warning: No connect value specified defaulting to \'False\'" False
    else
        let firstCommand = T.takeWhile (/= '=') $ head inTextList
            firstOption = T.tail $ T.dropWhile (/= '=') $ head inTextList
        in
        if firstCommand == T.pack "connect" then
            let option = T.unpack firstOption
            in
            if option == "true" then True
            else if option == "false" then False
            else error ("Connect option \'" ++ option ++ "\' not recognized (True|False)")
        else getConnect (tail inTextList)


-- | getThreshold returns threshold value or default otherwise
-- assumes in lower case
getThreshold :: [T.Text] -> Int
getThreshold inTextList =
    if null inTextList then trace "Warning: No threshold specified defaulting to \'0\'"  0 :: Int
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
    if null inTextList then trace "Warning: No output format specified defaulting to \'dot\'" "dot"
    else
        let firstCommand = T.takeWhile (/= '=') $ head inTextList
            firstOption = T.tail $ T.dropWhile (/= '=') $ head inTextList
        in
        if firstCommand == T.pack "outformat" then
            let outFormat = T.unpack firstOption
            in
            if outFormat == "dot" then "dot"
            else if outFormat == "fenewick" then "fennewick"
            else error ("Output format \'" ++ outFormat ++ "\' not recognized (dot|FENewickewick)")
        else getOutputFormat (tail inTextList)

-- | getOutputFileName returns output file name or default otherwise
-- assumes in lower case for command, uses pair so no case convewrsino in files name
getOutputFileName :: [(T.Text, T.Text)] -> String
getOutputFileName inTextPairList =
    if null inTextPairList then trace "Warning: No output file name specified defaulting to \'euncon.out\'" "euncon.out"
    else
        let (textListLC, textList) = head inTextPairList
            firstCommand = T.takeWhile (/= '=') textListLC
            firstOption = T.tail $ T.dropWhile (/= '=') textList
        in
        if firstCommand == T.pack "outfile" then T.unpack firstOption
        else getOutputFileName (tail inTextPairList)

