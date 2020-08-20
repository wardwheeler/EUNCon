{- |
Module      :  PhyloParsers.hs 
Description :  module witb parseing functios for commonly used phylogentic files
				graphs parsed to fgl types.
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

-}

module PhyloParse () where

import qualified Data.Graph.Inductive.Graph as G
import qualified Data.Graph.Inductive.PatriciaTree as P
import qualified Data.Text.Lazy as T


let newickTextList = filter (not.T.null) $ T.split (==';') $ removeNewickComments (T.concat newickFileTexts)
    -- hPutStrLn stderr $ show newickTextList
let newickGraphList = fmap (newickToGraph [] [] (-1, "") . (:[])) newickTextList


{--  
	Using Text as ouput for non-standard ascii charcaters (accents, umlautes etc)
--}


-- | getForestEnhancedNewickList takes String file contents and returns a list 
-- of fgl graphs with Text labels for nodes and edgesor error if not ForestEnhancedNewick or Newick formats.
getForestEnhancedNewickList :: String -> [P.Gr T.Text T.Text]
getForestEnhancedNewickList fleString = 
	if null fleString then error "Empty file string input in getForestEnhancedNewickList"
	else 
		let fileText = T.pack getForestEnhancedNewickList
			feNewickList = divideGraphText getForestEnhancedNewickList
		in
		fmap text2FGLGraph feNewickList

-- | divideGraphText splits multiple Text representations of graphs (Newick styles)
-- and returns a list if Text graph descriptions
-- also removed spaces from descriptions
-- converts 'Blah bleh" to Blah_bleh'
-- removes comments
divideGraphText :: T.Text -> [T.Text]
divideGraphText inText =
	if T.null inText then []
	else 
		[]

-- | removeBranchLengths from Text group
removeBranchLengths :: T.Text -> T.Text
removeBranchLengths inName
  | T.null inName = inName
  | T.last inName == ')' = inName
  | not (T.any (==':') inName) = inName
  | otherwise = T.reverse $ T.tail $ T.dropWhile (/=':') $ T.reverse inName
    -- )

-- | removeNewickComments take string and removes all "[...]"
removeNewickComments :: T.Text -> T.Text
removeNewickComments inString
  | T.null inString = T.empty
  | not (T.any (==']') inString) = inString
  | otherwise =
  let firstPart = T.takeWhile (/='[') inString
      secondPart = T.tail $ T.dropWhile (/=']') inString
  in
  T.append firstPart (removeNewickComments secondPart)

-- | text2FGLGraph takes Text of newick (forest or enhanced or OG) and
-- retns fgl graph representation
text2FGLGraph :: T.Text -> P.Gr T.Text T.Text
text2FGLGraph inGraphText = 
	if T.null inGraphText then error "Empty graph text in text2FGLGraph"
	else 
		let firstChar = T.head text2FGLGraph
			lastChar = T.last text2FGLGraph
		in
		if firstChar == '<' && lastChar == '>' then getFENewick inGraphText
		else if firstChar == '(' && lastChar == ';' then newickToGraph (T.init inGraphText)
		else error ("Graph text not in ForestEnhancedNewick or (Enhanced)Newick format")




-- | newickToGraph takes text of newick description (no terminal ';')
-- ther shoudl be no spaces in the text--filtered earlier.
-- and recurives creates a nodes and edges in fgl library style
-- progresively consumes the text, each apren block representing a node->subtree 
-- operates on lists of Text so can recurse over subgroups of whatever size
--NEED To KEEP NODE LABELS (especialky #) create labels for unlabelled nodes H-NUmber perhaps
newickToGraph :: [G.LNode T.Text] -> [G.LEdge T.Text] -> G.LNode T.Text -> [T.Text] -> P.Gr T.Text T.Text
newickToGraph nodeList edgeList parentNode inNewickTextList =
  -- trace ("Parsing " ++ (show inNewickTextList))  (
  if null inNewickTextList then
      -- remove edge to root and make graph
      G.mkGraph nodeList (filter ((> (-1)).fst3) edgeList)
  else
      let inNewickText = head inNewickTextList
      in
      if T.null inNewickText then G.empty
      -- is a leaf name only
      else if not (T.any (=='(') inNewickText) then
        let leaves = T.split (==',') inNewickText
            (newNodes, newEdges) = getLeafNodesAndEdges parentNode (length nodeList) leaves [] []
        in
        newickToGraph (nodeList ++ newNodes) (edgeList ++ newEdges) parentNode (tail inNewickTextList)
      -- should have '(' ')' and some charcters
      else if T.length inNewickText < 3 then error ("Improper newick tree component: " ++ T.unpack inNewickText)
      else
        --remove  labels (e.g. branch length as :XXX) on total group defined by '(...)' and remove outer parens
        let groupText = removeBranchLengths $ T.tail $ T.init $ removeBranchLengths inNewickText
            isSubTree = T.any (=='(') groupText
            thisNode = (length nodeList, "HTU" ++ show (length nodeList))
            thisEdge = (fst parentNode, length nodeList, "Edge(" ++ show (fst parentNode) ++ "," ++ show (length nodeList) ++ ")")
        in
        -- a set of leaves
        if not isSubTree then
          let leaves = T.split (==',') groupText
              (newNodes, newEdges) = getLeafNodesAndEdges thisNode (1 + length nodeList) leaves [] []
          in
          newickToGraph ((thisNode : nodeList) ++ newNodes) ((thisEdge : edgeList) ++ newEdges) parentNode (tail inNewickTextList)

        -- is a sub tree
        else
          -- find subtrees and recurse
          let subTrees = splitNewickGroups groupText
              -- spit into leaf and non-leaf descenedents
              (leafNodes, nonLeafNodes) = splitLeafNonLeafText subTrees [] []
              -- make any leaf nodes and edges
              (newNodes, newEdges) = getLeafNodesAndEdges thisNode (1 + length nodeList) leafNodes [] []
          in
          -- recurse the remainder
          newickToGraph ((thisNode : nodeList) ++ newNodes) ((thisEdge : edgeList) ++ newEdges) thisNode (nonLeafNodes ++ tail inNewickTextList)
          
