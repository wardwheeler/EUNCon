{- |
Module      :  GraphFormatUtilities.hs
Description :  module witb interconversion functions for commonly used phylogentic graph formats (newick. dot, fgl)
                graphs parsed to fgl types.
Copyright   :  (c) 2021 Ward C. Wheeler, Division of Invertebrate Zoology, AMNH. All rights reserved.
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


Forest Extended Newick defined here as a series of ENewick representations
within '<' ans '>'. Nodes can be shared among consituent ENewick representations
(';' from enewick itself, just for illustration, not doubled)

<EN1;EN2>

ExtendedNewick from Cardona et al. 2008.  BMC Bioinformatics 2008:9:532

    The labels and comments are as in Olsen Newick formalization below, except
    that underscores in unquoted label are NOT converted to spaces and quoted labels
    ar left as is with quotes and all.
    Other elements as in Cardona et all ENewick.



Gary Olsen's Interpretation of the "Newick's 8:45" Tree Format Standard
https://evolution.genetics.washington.edu/phylip/newick_doc.html

Conventions:
   Items in { } may appear zero or more times.
   Items in [ ] are optional, they may appear once or not at all.
   All other punctuation marks (colon, semicolon, parentheses, comma and
         single quote) are required parts of the format.

              tree ==> descendant_list [ root_label ] [ : branch_length ] ;

   descendant_list ==> ( subtree { , subtree } )

           subtree ==> descendant_list [internal_node_label] [: branch_length]
                   ==> leaf_label [: branch_length]

            root_label ==> label
   internal_node_label ==> label
            leaf_label ==> label

                 label ==> unquoted_label
                       ==> quoted_label

        unquoted_label ==> string_of_printing_characters
          quoted_label ==> ' string_of_printing_characters '

         branch_length ==> signed_number
                       ==> unsigned_number

Notes:
   Unquoted labels may not contain blanks, parentheses, square brackets,
        single_quotes, colons, semicolons, or commas.
   Underscore characters in unquoted labels are converted to blanks.
   Single quote characters in a quoted label are represented by two single
        quotes.
   Blanks or tabs may appear anywhere except within unquoted labels or
        branch_lengths.
   Newlines may appear anywhere except within labels or branch_lengths.
   Comments are enclosed in square brackets and may appear anywhere
        newlines are permitted.

Other notes:
   PAUP (David Swofford) allows nesting of comments.
   TreeAlign (Jotun Hein) writes a root node branch length (with a value of
        0.0).
   PHYLIP (Joseph Felsenstein) requires that an unrooted tree begin with a
        trifurcation; it will not "uproot" a rooted tree.

Example:
   (((One:0.2,Two:0.3):0.3,(Three:0.5,Four:0.3):0.2):0.3,Five:0.7):0.0;

           +-+ One
        +--+
        |  +--+ Two
     +--+
     |  | +----+ Three
     |  +-+
     |    +--+ Four
     +
     +------+ Five
--}

{- 
ToDo: 
    Need to verify FEN output and input with Dendroscope (errors/non-parse reported)
-}

module GraphFormatUtilities (forestEnhancedNewickStringList2FGLList,
                             fglList2ForestEnhancedNewickString,
                             component2Newick,
                             checkIfLeaf,
                             stringGraph2TextGraph,
                             textGraph2StringGraph,
                             showGraph,
                             relabelFGL,
                             convertGraphToStrictText,
                             splitVertexList,
                             relabelFGLEdgesDouble
                            ) where

import           Control.Parallel.Strategies
import           Data.Char                         (isSpace)
import qualified Data.Graph.Inductive.Graph        as G
import qualified Data.Graph.Inductive.PatriciaTree as P
import           Data.List
import qualified Data.Map.Strict                   as Map
import           Data.Maybe
import qualified Data.Text.Lazy                         as T
import qualified Data.Text                         as StrictT

import           Debug.Trace
import           GeneralUtilities
import           ParallelUtilities
import           Data.GraphViz                     as GV
import           Data.GraphViz.Attributes.Complete (Attribute (Label),
                                                    Label (..))
import           Data.GraphViz.Commands.IO         as GVIO
import           Data.Monoid



-- | showGraph a semi-formatted show for Graphs
showGraph :: (Show a, Show b) => P.Gr a b -> String -- BV.BV (BV.BV, BV.BV) -> String
showGraph inGraph =
  if G.isEmpty inGraph then "Empty Graph"
  else
      let nodeString = show $ G.labNodes inGraph
          edgeString  = show $ G.labEdges inGraph
      in
      ("Nodes:" ++ nodeString ++ "\n" ++ "Edges: " ++ edgeString)

-- | getForestEnhancedNewickList takes String file contents and returns a list
-- of fgl graphs with Text labels for nodes and edges or error if not ForestEnhancedNewick or Newick formats.
forestEnhancedNewickStringList2FGLList :: T.Text -> [P.Gr T.Text Double]
forestEnhancedNewickStringList2FGLList fileText =
    if T.null fileText then []
    else
        let feNewickList = fmap (removeNewickSpaces . removeNewickComments) (divideGraphText fileText)
        in
        -- trace ("There are " ++ (show $ length feNewickList) ++ " graphs to convert: " ++ (show feNewickList))
        parmap rdeepseq text2FGLGraph feNewickList

-- | divideGraphText splits multiple Text representations of graphs (Newick styles)
-- and returns a list of Text graph descriptions
divideGraphText :: T.Text -> [T.Text]
divideGraphText inText =
    if T.null inText then []
    else
        let firstChar = T.head inText
        in
        if firstChar == '<' then
            let firstPart = T.snoc (T.takeWhile (/= '>') inText) '>'
                restPart = T.tail $ T.dropWhile (/= '>') inText
            in
            firstPart : divideGraphText restPart
        else if firstChar == '(' then
            let firstPart = T.snoc (T.takeWhile (/= ';') inText) ';'
                restPart = T.tail $ T.dropWhile (/= ';') inText
            in
            firstPart : divideGraphText restPart
        else error ("First character in graph representation " ++ (T.unpack inText) ++ " : " ++ (show firstChar) ++ " is not either '<' or '('")

-- | removeNewickComments take text and removes all "[...]"
removeNewickComments :: T.Text -> T.Text
removeNewickComments inString
  | T.null inString = T.empty
  | not (T.any (==']') inString) = inString
  | otherwise =
  let firstPart = T.takeWhile (/='[') inString
      secondPart = T.tail $ T.dropWhile (/=']') inString
  in
  T.append firstPart (removeNewickComments secondPart)

-- | convertQuotedText takes single quoted Text and removes single quotes and converts spaces
-- to underscores
convertQuotedText :: T.Text -> (T.Text, T.Text)
convertQuotedText inText =
  if T.null inText then error "Emmpty Text in convertQuotedText"
  else
    let firstPart = T.replace (T.singleton ' ') (T.singleton '_') $ T.takeWhile (/= '\'') (T.tail inText)
        restPart = T.tail $ T.dropWhile (/= '\'') (T.tail inText)
    in
    (firstPart, restPart)

-- | removeNewickSpaces removes spaces and converts single quoted strings
-- with spaces to unquoted strings with underscores replacing spaces: ' blah bleh ' => blah_bleh
removeNewickSpaces :: T.Text -> T.Text
removeNewickSpaces inText =
  if T.null inText then T.empty
  else
    let firstChar = T.head inText
    in
    if firstChar == '\'' then
      let (newText, restText) = convertQuotedText inText
      in
      T.concat [newText, removeNewickSpaces restText]
    else if isSpace firstChar then removeNewickSpaces $ T.tail inText
    else T.cons firstChar (removeNewickSpaces $ T.tail inText)

-- | text2FGLGraph takes Text of newick (forest or enhanced or OG) and
-- retns fgl graph representation
text2FGLGraph :: T.Text -> P.Gr T.Text Double
text2FGLGraph inGraphText =
    if T.null inGraphText then error "Empty graph text in text2FGLGraph"
    else
        let firstChar = T.head inGraphText
            lastChar = T.last inGraphText
        in
        if firstChar == '<' && lastChar == '>' then fENewick2FGL inGraphText -- getFENewick inGraphText
        else if firstChar == '(' && lastChar == ';' then mergeNetNodesAndEdges $ makeGraphFromPairList $ eNewick2FGL [] []  [(inGraphText, (-1, T.empty))]
        else error "Graph text not in ForestEnhancedNewick or (Enhanced)Newick format"


-- | fENewick2FGL takes a Forest Extended Newick (Text) string and returns FGL graph
-- breaks up forest and parses seprate eNewicks then modifes for any
-- common network nodes in the sub-graphs
fENewick2FGL :: T.Text -> P.Gr T.Text Double
fENewick2FGL inText =
  if T.null inText then error "Empty graph text in fENewick2FGL"
    else
      -- split eNewicks
      let eNewickTextList = splitForest inText
          startNodeList = replicate (length eNewickTextList) (-1, T.empty)
          textNodeList = zip eNewickTextList startNodeList
      -- init to remove trailing ';' from eNewick
          eNewickGraphList = fmap ((mergeNetNodesAndEdges . makeGraphFromPairList) . (eNewick2FGL [] [] . (:[]))) textNodeList
      in
      if length eNewickGraphList == 1 then head eNewickGraphList
      else
          -- merge graphs then merge network nodes and edges edges
          let fENewickGraph = mergeNetNodesAndEdges $ mergeFGLGraphs G.empty eNewickGraphList
          in
          fENewickGraph

-- | splitForest takes a Text (string) Forest Enhanced Newick representation and splits into
-- its consituent Extended Newick representations
splitForest :: T.Text -> [T.Text]
splitForest inText
  | T.null inText = []
  | (T.head inText /= '<') || (T.last inText /= '>') = error ("Invalid Forest Extended Newick representation," ++
    " must begin with \'<\'' and end with \'>\' : " ++ T.unpack inText)
  | otherwise =
  let partsList = filter (not.T.null) $ T.splitOn (T.singleton ';') (T.init $ T.tail inText)
      eNewickList = fmap (`T.append` T.singleton ';') partsList
  in
  eNewickList

-- | makeGraphFromPairList takes pair of node list and edge list and returns Graph
-- | filters to remove place holder node and edges creted during eNewick pass
makeGraphFromPairList :: [(G.LNode T.Text,G.LEdge Double)] -> P.Gr T.Text Double
makeGraphFromPairList pairList =
  if null pairList then G.empty
  else
    let (nodeList, edgeList) = unzip pairList
    in
    G.mkGraph (filter ((> (-1)).fst) nodeList) (filter ((> (-1)).fst3) edgeList)

-- | getBranchLength extracts branch length from Text label and puts in '1' if there is no
-- branch length--makes sure after last ')'
getBranchLength :: T.Text -> Double
getBranchLength inText =
  --trace ("Getting branch length of " ++ show inText) (
  if T.null inText then error "Null text in getBranchLength"
  else
    let a = T.dropWhile (/= ':') $ T.reverse $ T.takeWhile (/= ')') $ T.reverse inText
    in
    if T.null a then 1
    else if T.length a == 1 then error "Need branch length after \':\')"
    else (read (T.unpack $ T.tail a) :: Double)
    --)

-- | getNodeLabel get--or makes--a label for a node
-- after last ')' before any ':', without ',' after last ')'
getNodeLabel :: Int -> T.Text -> T.Text
getNodeLabel nodeNumber inText =
  --trace ("Getting node label of " ++ show inText) (
  if T.null inText then error "Null text in getNodeLabel"
  else
    let a = T.takeWhile (/= ':') $ T.reverse $ T.takeWhile (/= ')') $ T.reverse inText
    in
    if T.any (==',') a || T.null a then T.append (T.pack "HTU") (T.pack $ show nodeNumber) else a
    --)

-- | getLeafInfo takes Text of teminal (no ',') and parses to yeild
-- either a single leaf label, edge, and edge weight, or two
-- leaves with labels and costs if there is a network node as parent
-- need to merge network nodes later
getLeafInfo :: T.Text -> G.LNode T.Text -> [G.LNode T.Text] -> [(G.LNode T.Text,G.LEdge Double)]
getLeafInfo leafText parentNode nodeList
  | T.null leafText = error "Empty leaf text in getLeafInfo"
  | not (T.any (=='(') leafText) =
  let leafLabel = T.takeWhile (/= ':') leafText
      edgeWeight = getBranchLength leafText
      -- CHECK FOR EXISTING
      thisNode = (length nodeList, leafLabel)
      thisEdge = (fst parentNode, length nodeList, edgeWeight)
      --preexistingNode = checkForExistingNode leafLabel nodeList
  in
  [(thisNode, thisEdge)]
  | otherwise =
  let -- leaf parent info
      -- (leafLabel)leafParentLabel:leafParentBranchLength
      leafParentEdgeWeight = getBranchLength leafText
      leafParentLabel = getNodeLabel (length nodeList) leafText
      leafParentNode = (length nodeList, leafParentLabel)
      leafParentEdge = (fst parentNode, fst leafParentNode, leafParentEdgeWeight)

      -- leaf info
        -- (leafLabel)X#H:000 => leafLabel
      leafLabelText = T.takeWhile (/= ')') $ T.tail leafText
      -- check for existing
      leafLabel = T.takeWhile (/= ':') leafLabelText
      leafEdgeWeight = getBranchLength leafLabelText
      leafNode = (1 + length nodeList, leafLabel)
      leafEdge = (fst leafParentNode, fst leafNode, leafEdgeWeight)
  in
  [(leafNode, leafEdge),(leafParentNode, leafParentEdge)]

-- | getBodyParts takes a Text of a subTree and splits out the group description '(blah)', any node label
-- and any branch length
getBodyParts :: T.Text -> Int -> (T.Text, T.Text, Double)
getBodyParts inRep nodeNumber =
  if T.null inRep then error "No group to parse in getBodyParts"
  else
      --trace ("In body parts") (
      let subGraphPart =  T.reverse $ T.dropWhile (/= ')') $ T.reverse inRep
          branchLength =  getBranchLength inRep
          subGraphLabel = getNodeLabel nodeNumber inRep
      in
      --trace (show (subGraphPart, subGraphLabel, branchLength))
      (subGraphPart, subGraphLabel, branchLength)
      --)

-- | getParenBoundedGraph tkaes a Text String and returns  the first graph component
-- with balanced parens and remainder of Text
getParenBoundedGraph :: Int -> Int -> T.Text -> T.Text -> (T.Text, T.Text)
getParenBoundedGraph leftParenCounter rightParenCounter curText inText =
  --trace ("GB " ++ show curText ++ " " ++ show inText) (
  if T.null inText then (curText, inText)
  else
    let firstChar = T.head inText
    in
    if firstChar == '(' then getParenBoundedGraph (leftParenCounter + 1) rightParenCounter (T.snoc curText firstChar) (T.tail inText)
    else if firstChar /= ')' then getParenBoundedGraph leftParenCounter rightParenCounter (T.snoc curText firstChar) (T.tail inText)
    else -- right paren
      if rightParenCounter + 1 == leftParenCounter then -- closing matrched paren
          let restOfComponent = T.takeWhile (/= ',') inText
              remainderText = T.dropWhile (/= ',') inText
          in
          (curText `T.append` restOfComponent, remainderText)
      else getParenBoundedGraph leftParenCounter (rightParenCounter + 1) (T.snoc curText firstChar) (T.tail inText)
     -- )

-- | getSubComponents takes a Text String and reurns 1 or more subcomponents of graph
-- scenarios include leaf, leaf in parens, subgraph in parens
getSubComponents :: T.Text -> [T.Text]
getSubComponents inText
  | T.null inText = []
  | T.head inText == ',' = getSubComponents (T.tail inText)
  | T.head inText /= '(' = -- simple leaf (no net node labels)
  let subGraph = T.takeWhile (/= ',') inText
      restGraph = T.dropWhile (/= ',') inText
  in
  subGraph : getSubComponents restGraph
  | otherwise = -- "regular" paren defined element
  let (subGraph, restGraph) = getParenBoundedGraph 0 0 T.empty inText
  in
  subGraph : getSubComponents restGraph

-- | getChildren splits a subGraph Text '(blah, blah)' by commas, removing outer parens
getChildren :: T.Text -> [T.Text]
getChildren inText
  | T.null inText = []
  | (T.head inText /= '(') || (T.last inText /= ')') = error ("Invalid Extended Newick component," ++
    " must begin with \'(\'' and end with \')\' : " ++ T.unpack inText)
  | otherwise =
  -- modify for indegree 1 outdegree 1 and print warning.
  let guts = T.init $ T.tail inText -- removes leading and training parens
      subComponents = filter (not.T.null) $  getSubComponents guts
  in
  subComponents

-- | checkForExistingNode takes a node label and checs the node list for the first
-- node with the same label and returns a Maybe node, else Nothing
checkForExistingNode :: T.Text -> [G.LNode T.Text] -> Maybe (G.LNode T.Text)
checkForExistingNode nodeLabel nodeList =
  if null nodeList then Nothing
  else
    let matchList = filter ((==nodeLabel).snd) nodeList
    in
    if null matchList then Nothing
    else Just $ head matchList

-- | checkIfLeaf checks text to see if leaf.
-- if the number of left parens is 1 and right parens 1 and no ',' then leaf
-- if no left parens and no right paren then leaf
-- then its a leaf
-- either "bleh", "bleh:00", or "(bleh)label:00"
checkIfLeaf :: T.Text -> Bool
checkIfLeaf inText =
  if T.null inText then error "Null text to check if leaf in checkIfLeaf"
  else
    let leftParenCount = T.count (T.singleton '(') inText
        rightParenCount = T.count (T.singleton ')') inText
        commaCount = T.count (T.singleton ',') inText
    in
    ((leftParenCount == 0) && (rightParenCount == 0) && (commaCount == 0)) || (if (leftParenCount == 0) && (rightParenCount == 0) && (commaCount > 0) then error ("Comma within leaf label" ++ show inText)
                                                                          else (leftParenCount == 1) && (rightParenCount == 1) && (commaCount == 0))

-- | eNewick2FGL takes a single Extended Newick (Text) string and returns FGL graph
-- allows arbitrary in and out degree except for root and leaves
eNewick2FGL :: [G.LNode T.Text] -> [G.LEdge Double] -> [(T.Text, G.LNode T.Text)] -> [(G.LNode T.Text,G.LEdge Double)]
eNewick2FGL nodeList edgeList inTextParentList =
    if null inTextParentList then []
    else
      let inTextFirst = fst $ head inTextParentList
          parentNode = snd $ head inTextParentList
          isRoot = null nodeList
      in
      -- see if initial call and check format
      if isRoot && ((T.head inTextFirst /= '(') || (T.last inTextFirst /= ';'))  then error ("Invalid Extended Newick component," ++
      " must begin with \'(\'' and end with \')\' : " ++ T.unpack inTextFirst)
      -- not first call and/or format OK
      else
        let inText = if isRoot then T.takeWhile (/= ';') inTextFirst else inTextFirst -- remove trailing ';' if first (a bit wasteful--but intial check on format)
            isLeaf = checkIfLeaf inText
        in
        -- trace ("Parsing " ++ show inText ++ " from parent " ++ show parentNode ++ " " ++ show isLeaf)(
        -- is a single leaf
          -- need better could be  series of indegree `1 outdegree 1 nodes to a single leaf with no ','
        -- ike (a(b(c(d))))
        if isLeaf then
          -- parse label ala Gary Olsen formalization
          -- since could have reticulate label yeilding two edges and two nodes
          -- Cardona et al 2008  Extended Newick
          let newLeafList = getLeafInfo inText parentNode nodeList
              newNodeList = fmap fst newLeafList
              newEdgeList = fmap snd newLeafList
          in
          newLeafList ++ eNewick2FGL (newNodeList ++ nodeList) (newEdgeList ++ edgeList) (tail inTextParentList)
        else
          -- is subtree assumes start and end with parens '(blah)'
          let (subTree, nodeLabel, edgeWeight) = getBodyParts inText (length nodeList)
              thisNode = (length nodeList, nodeLabel)
              thisEdge = (fst parentNode, length nodeList, edgeWeight)
              childTextList = getChildren subTree
              parentNodeList = replicate (length childTextList) thisNode
              childParentList = zip childTextList parentNodeList

          in
          (thisNode, thisEdge) : eNewick2FGL (thisNode : nodeList) (thisEdge : edgeList) (childParentList ++ tail inTextParentList)

-- | reindexNode takes an offset and adds to the node index
-- returning new node
reindexNode :: Int -> G.LNode T.Text -> G.LNode T.Text
reindexNode offSet (index, label) = (index + offSet, label)

-- | reindexEdge takes an offset and adds to the two indices of the edge
-- returning the new edge
reindexEdge :: Int -> G.LEdge Double -> G.LEdge Double
reindexEdge offSet (e, u, label) = (e + offSet, u + offSet, label)

-- | mergeFGLGraphs takes multiple graphs and merges
-- nodes and edges via reindexing
-- just adds progessive offsets from graph node indices as added
mergeFGLGraphs :: P.Gr T.Text Double -> [P.Gr T.Text Double] -> P.Gr T.Text Double
mergeFGLGraphs curGraph inGraphList
  | null inGraphList = curGraph
  | G.isEmpty curGraph = mergeFGLGraphs (head inGraphList) (tail inGraphList)
  | otherwise =
  let firstGraph = head inGraphList
      firstNodes = G.labNodes firstGraph
      firstEdges = G.labEdges firstGraph
      curNodes = G.labNodes curGraph
      curEdges = G.labEdges curGraph
      newNodes = fmap (reindexNode (length curNodes)) firstNodes
      newEdges = fmap (reindexEdge (length curNodes)) firstEdges
  in
  mergeFGLGraphs (G.mkGraph (curNodes ++ newNodes) (curEdges ++ newEdges)) (tail inGraphList)

-- | getNodeIndexPair take a list of unique nodes and checks successive nodes and
-- adds to unique list, also creating a full list of pairs of indicess for non-unique that
-- can be used as an index map for edges
-- length of unique list to keep the node indices sequential
getNodeIndexPair :: [G.LNode T.Text] -> [(Int, Int)] -> [G.LNode T.Text] -> ([G.LNode T.Text], [(Int, Int)])
getNodeIndexPair uniqueList pairList nodeToCheckList =
  if null nodeToCheckList then (reverse uniqueList, reverse pairList)
  else
    let firstNode@(index, label) = head nodeToCheckList
        matchingNode = checkForExistingNode label uniqueList
    in
    if isNothing matchingNode then getNodeIndexPair (firstNode : uniqueList) ((index, length uniqueList) : pairList) (tail nodeToCheckList)
    else
      let existingNode = fromJust matchingNode
          newPair = (index, fst existingNode)
      in
      getNodeIndexPair uniqueList (newPair : pairList) (tail nodeToCheckList)

-- | reIndexEdge takes an (Int, Int) map, labelled edge, and returns a new labelled edge with new e,u vertices
reIndexLEdge ::  Map.Map Int Int -> G.LEdge Double -> G.LEdge Double
reIndexLEdge vertexMap inEdge =
  if Map.null vertexMap then error "Null vertex map"
  else
    let (e,u,label) = inEdge
        newE = Map.lookup e vertexMap
        newU = Map.lookup u vertexMap
    in
    if isNothing newE then error ("Error looking up vertex " ++ show e ++ " in " ++ show (e,u))
    else if isNothing newU then error ("Error looking up vertex " ++ show u ++ " in " ++ show (e,u))
    else (fromJust newE, fromJust newU, label)

-- | reIndexNode takes an (Int, Int) map, labelled node, and returns a new labelled node with new vertex
reIndexLNode ::  Map.Map Int Int -> G.LNode T.Text -> G.LNode T.Text
reIndexLNode vertexMap inNode =
  if Map.null vertexMap then error "Null vertex map"
  else
    let (index,label) = inNode
        newIndex = Map.lookup index vertexMap
    in
    if isNothing newIndex then error ("Error looking up vertex " ++ show index ++ " in " ++ show inNode)
    else (fromJust newIndex, label)

-- | mergeNetNodesAndEdges takes a single graph and merges
-- nodes and edges due to network nodes and edges
-- uses checkForExistingNode and creates a map from nodes to reindex edges
-- needs to be merged first if graphs are combined--or indices will be wrong
mergeNetNodesAndEdges :: P.Gr T.Text Double -> P.Gr T.Text Double
mergeNetNodesAndEdges inGraph =
  if G.isEmpty inGraph then G.empty
  else
    let nodeList = G.labNodes inGraph
        (_, nodeIndexPairs) = getNodeIndexPair [] [] nodeList
        nodeMap = Map.fromList nodeIndexPairs
        reindexedNodeList = fmap (reIndexLNode nodeMap) (G.labNodes inGraph)
        reIndexedEdgeList = fmap (reIndexLEdge nodeMap) (G.labEdges inGraph)
    in
    G.mkGraph reindexedNodeList reIndexedEdgeList

-- | subTreeSize takes a nodeList and retuns the number of leaves that can be
-- traced back to those nodes (for single just pass list of 1 node)
-- this used for ordering of groups left (smaller) right (larger)
subTreeSize :: P.Gr a b -> Int -> [G.LNode a] -> Int
subTreeSize inGraph counter nodeList =
  if null nodeList then counter
  else
    let firstNode = head nodeList
        children = G.suc inGraph $ fst firstNode
        -- assumes all childrfen have a label (if not--problems)
        labelList = fmap (fromJust . G.lab inGraph) children
        labChildren = zip children labelList
    in
    subTreeSize inGraph (counter + length labChildren) (tail nodeList ++ labChildren)

-- | getRoot takes a greaph and list of nodes and returns vertex with indegree 0
-- so assumes a connected graph--with a single root--not a forest
getRoots :: P.Gr a b -> [G.LNode a] -> [G.LNode a]
getRoots inGraph nodeList =
  if null nodeList then [] --error "Root vertex not found in getRoot"
  else
    let firstNode@(index, _) = head nodeList
    in
    if G.indeg inGraph index == 0 then firstNode : getRoots inGraph (tail nodeList)
    else getRoots inGraph (tail nodeList)

-- | removeDuplicateSubtreeText removes duplicate subtree textx that come from indegree > 1 nodes
-- there should be at least two of each network texts.
-- for each case, the first instance is kept, and the remainders are replaced with the node label
-- and edge weight if specified (:000)
removeDuplicateSubtreeText :: (Show b) => T.Text -> [G.LNode T.Text] -> P.Gr T.Text b -> Bool -> Bool -> T.Text
removeDuplicateSubtreeText inRep netNodeList fglGraph writeEdgeWeight writeNodeLable =
  if null netNodeList then inRep
  else
    let netNodeText = T.init $ component2Newick fglGraph writeEdgeWeight writeNodeLable (head netNodeList)
        -- edge weight already removed or may not match all occurences
        -- I have no idea why--but there are extraneous double quotes that have to be removed.
        nodeText' = T.filter (/= '\"') netNodeText
        -- checks to see if leaf is indegree > 1
        nodeText = if not (checkIfLeaf nodeText') then nodeText' else T.filter (/= '(') $ T.filter (/= ')') nodeText'
        nodeLabel = T.reverse $ T.takeWhile (/= ')') $ T.reverse nodeText
        textList = T.splitOn nodeText inRep
        -- isFound = T.isInfixOf nodeText inRep -- (T.pack "(4:1.0)Y#H1") inRep
    in
    -- trace ("Removing ? " ++ show isFound ++ " " ++ show nodeText ++ " " ++ show nodeLabel ++ " from " ++ show inRep ++ " in list (" ++ show (length textList) ++ ") " ++ show textList) (
    -- since root cannot be network neither first nor last pieces should be empty
    if T.null (head textList) || T.null (last textList) then error ("Text representation of graph is incorrect with subtree:\n" ++ T.unpack nodeText
      ++ " first or last in representation:\n " ++ T.unpack inRep)
    --else if length textList == 1 then error ("Text representation of graph is incorrect with subtree:\n" ++ T.unpack nodeText
    --  ++ " not found in representation:\n " ++ T.unpack inRep)
    else if length textList == 2 then
        trace "Warning: Network subtree present only once--extraneous use of \'#\' perhaps--"
        inRep
    else -- should be minimum of 3 pieces (two occurences of subtree text removed) is subtree found 2x and not sister to itself (which should never happen)
      -- edge weights (if they occur) remain the beginning of each text list (after the first)
      let firstPart = head textList `T.append` nodeText
          secondPart = T.intercalate nodeLabel (tail textList)
      in
      removeDuplicateSubtreeText (firstPart `T.append` secondPart) (tail netNodeList) fglGraph writeEdgeWeight writeNodeLable
      --)

-- | getDistToRoot takes a node and a graph and gets the shortes path to root
-- and returns the number of links
getDistToRoot :: P.Gr T.Text b -> Int -> G.Node -> Int
getDistToRoot fglGraph counter inNode  =
  if counter > length (G.nodes fglGraph) then error "Cycle likely in graph, path to root larger than number of nodes"
  else
    let parents = G.pre fglGraph inNode
    in
    if null parents then counter
    else
      let parentPaths = fmap (getDistToRoot fglGraph (counter + 1)) parents
      in
      minimum parentPaths

-- | modifyInDegGT1Leaves operates on the leaf nodes with indegree > 1 to prepare them for enewick representation
-- leaves with indegree greater than one are modified such that:
--   1) a new node is created as parent to the leaf and added to the "add" list
--   2) edges incident on the leaf are put in the "delete" list 
--   3) edges are created from teh new node to the leaf, and edges are "added" from teh leaf's parent to the new node
modifyInDegGT1Leaves :: P.Gr T.Text Double -> Int -> [G.LNode T.Text] -> ([G.LNode T.Text],[G.LEdge Double],[G.LEdge Double]) -> ([G.LNode T.Text],[G.LEdge Double],[G.LEdge Double]) 
modifyInDegGT1Leaves origGraph totalNumberNodes leafNodes graphDelta@(nodesToAdd, edgesToAdd, edgesToDelete) =
  if G.isEmpty origGraph then error "Empty graph in modifyInDegGT1Leaves"
  else if null leafNodes then graphDelta
  else 
      let firstLeaf = head leafNodes
          firstInDegree = G.indeg origGraph (fst firstLeaf)
      in
      if firstInDegree == 1 then modifyInDegGT1Leaves origGraph totalNumberNodes (tail leafNodes) graphDelta
      else -- in a "network leaf"
          let inEdgeList = G.inn origGraph (fst firstLeaf)
              parentNodeList = fmap fst3 inEdgeList
              inEdgeLabels = fmap thd3 inEdgeList
              newNode = (totalNumberNodes, T.pack $ "HTU" ++ show totalNumberNodes)
              repeatedNodeNumber = replicate (length inEdgeList) totalNumberNodes
              newEdgeList = (totalNumberNodes, (fst firstLeaf), 0 :: Double) : (zip3 parentNodeList repeatedNodeNumber inEdgeLabels)
          in
          modifyInDegGT1Leaves origGraph (totalNumberNodes + 1) (tail leafNodes) (newNode : nodesToAdd, newEdgeList, inEdgeList)

-- | modifyInDegGT1HTU operates on the HTU nodes with indegree > 1 to prepare them for enewick representation
-- HTUs with indegree greater than one are modified such that:
--   1) the original HTU is maintained the first edge to that HTU is Maintained also
--   2) other edges to that HTU are put in "delete" list
--   3) a new HTU node is created for each deleted edge with same label as original HTU
--   4) edges are created from the parents (except for the firt one) to the new nodes such thayt each one is indegree=1 and outdegree=0
--   5) ne edges are put in teh "add list"
modifyInDegGT1HTU :: P.Gr T.Text Double -> Int -> [G.LNode T.Text] -> ([G.LNode T.Text],[G.LEdge Double],[G.LNode T.Text],[G.LEdge Double]) -> ([G.LNode T.Text],[G.LEdge Double],[G.LNode T.Text],[G.LEdge Double]) 
modifyInDegGT1HTU origGraph totalNumberNodes htuNodes graphDelta@(nodesToAdd, edgesToAdd, nodesToDelete, edgesToDelete) =
  if G.isEmpty origGraph then error "Empty graph in modifyInDegGT1HTU"
  else if null htuNodes then graphDelta
  else 
      let firstHTU = head htuNodes
          firstInDegree = G.indeg origGraph (fst firstHTU)
      in
      if firstInDegree == 1 then modifyInDegGT1HTU origGraph totalNumberNodes (tail htuNodes) graphDelta
      else -- in a "network leaf"
          let inEdgeList = G.inn origGraph (fst firstHTU)
              parentNodeList = fmap fst3 inEdgeList
              inEdgeLabels = fmap thd3 inEdgeList

              -- Create new nodes and edges
              newNodeList = []
              newEdgeList = edgesToAdd
              {-
              newNodeList = (totalNumberNodes, T.pack $ "HTU" ++ show totalNumberNodes)
              repeatedNodeNumber = replicate (length inEdgeList) totalNumberNodes
              newEdgeList = (totalNumberNodes, (fst firstLeaf), 0 :: Double) : (zip3 parentNodeList repeatedNodeNumber inEdgeLabels)
              -}
          in
          modifyInDegGT1HTU origGraph (totalNumberNodes + 1) (tail htuNodes) (newNodeList ++ nodesToAdd, newEdgeList, nodesToDelete, inEdgeList)

-- | modifyFGLForEnewick takes an FGl graphs and modified for enewick output
-- 1) makes leaves that are indegree > 1 indegree 1 by creationg a new parent node with 
--    the leafs parents as parents and teh leaf as single child
-- 2) convertes each indegree > 1 node (non--leaf) to a series of indegree 1 nodes
--    the first of which gets the first parent and all the children. Each subsequent
--    parent gets a node with the name label and no children
modifyFGLForEnewick :: P.Gr T.Text Double -> P.Gr T.Text Double
modifyFGLForEnewick inGraph =
  if G.isEmpty inGraph then error "Empty graph to modify in modifyFGLForEnewick"
  else
    let inNodes = G.labNodes inGraph
        inEdges = G.labEdges inGraph
        (rootNodes, leafNodes, nonLeafNodes) = splitVertexList inGraph
        (nodesToAdd, edgesToAdd, edgesToDelete) = modifyInDegGT1Leaves inGraph (length inNodes) leafNodes ([],[],[])
        leafModGraph = G.insEdges edgesToAdd $ G.insNodes nodesToAdd $ G.delEdges (fmap G.toEdge edgesToDelete) inGraph
        (nodesToAddHTU, edgesToAddHTU, nodesToDelete, edgesToDeleteHTU) = modifyInDegGT1HTU leafModGraph (G.order leafModGraph) nonLeafNodes ([],[],[],[])
    in
    G.insEdges edgesToAddHTU $ G.insNodes nodesToAddHTU $ G.delNodes (fmap fst nodesToDelete) $ G.delEdges (fmap G.toEdge edgesToDeleteHTU) leafModGraph
        

-- | fgl2FEN take a fgl graph and returns a Forest Enhanced Newick Text
--   Can be simplified along lines of Cardona with change to graph before 
--   generating the rep bu splitting Hybrid nodes.  
--   enewick requires (it seems) indegree 1 for leaves
--   so need to creates nodes for indegree > 1 leaves
-- these are not issues for dot files
fgl2FEN :: Bool -> Bool -> P.Gr T.Text Double -> T.Text
fgl2FEN writeEdgeWeight writeNodeLable inFGLGraph =
  if G.isEmpty inFGLGraph then error "Empty graph to convert in fgl2FEN"
  else
    -- Modify greaph for enewick stuff (leaves -> indegree 1, 'split' network nodes)
    let fglGraph = modifyFGLForEnewick inFGLGraph
    in
    -- get forest roots
    let numRoots = getRoots fglGraph (G.labNodes fglGraph)
        rootGraphSizeList = fmap (subTreeSize fglGraph 0 . (:[])) numRoots
        rootAndSizes = zip rootGraphSizeList numRoots
        rootOrder = sortOn fst rootAndSizes
        fenTextList = fmap (component2Newick fglGraph writeEdgeWeight writeNodeLable . snd) rootOrder
        wholeRep = T.concat $ (`T.append` T.singleton '\n') <$> fenTextList

        -- filter stuff for repeated subtrees
        networkNodeList = filter ((>1).G.indeg fglGraph) (G.nodes fglGraph) -- not labNodes becasue of indeg

        -- need to order based on propinquity to root so fixes nested network issues
        -- closer to root first
        distToRootList = fmap (getDistToRoot fglGraph 0) networkNodeList
        netNodeRootPairsList = zip distToRootList networkNodeList
        orderedNetworkNodeList = snd <$> sortOn fst netNodeRootPairsList

        labelList = fmap (makeLabel . G.lab fglGraph) orderedNetworkNodeList
        networkLabelledNodeList = zip orderedNetworkNodeList labelList
        {-
        Moving netNodeTreeTextList into remove duplicate because can be nested.
        -- remove trailing ';''
        netNodeTreeTextList = fmap T.init $ fmap (component2Newick fglGraph writeEdgeWeight writeNodeLable) networkLabelledNodeList
        wholeRep' = removeDuplicateSubtreeText wholeRep netNodeTreeTextList
        -}
        -- this can be removed if graph modified as in Cardona
        wholeRep' = removeDuplicateSubtreeText wholeRep networkLabelledNodeList fglGraph writeEdgeWeight writeNodeLable
    in
    -- trace ("fgl2FEN " ++ (show $ length numRoots) ++ " " ++ show rootOrder ++ "->" ++ show fenTextList) (
    if length fenTextList == 1 then wholeRep' -- just a single tree/network
    else T.snoc (T.cons '<' wholeRep') '>' -- forest
    -- )


-- | fglList2ForestEnhancedNewickString takes FGL representation of forest and returns
-- list of Forest Enhanced Newick as a single String
fglList2ForestEnhancedNewickString :: [P.Gr T.Text Double] -> Bool -> Bool -> String
fglList2ForestEnhancedNewickString inFGLList writeEdgeWeight writeNodeLable =
  if null inFGLList then "\n"
  else
    let forestTextList = (`T.append` T.singleton '\n') <$> parmap rdeepseq (fgl2FEN writeEdgeWeight writeNodeLable) inFGLList
        forestListString = T.unpack $ T.concat forestTextList
    in
    forestListString

-- | component2Newick take a graph and root and creates enhanced newick from that root
component2Newick :: (Show a) => P.Gr T.Text a -> Bool -> Bool -> G.LNode T.Text -> T.Text
component2Newick fglGraph writeEdgeWeight writeNodeLable (index, label) =
  if G.isEmpty fglGraph then error "Empty graph to convert in component2Newick"
  else
    -- start with root (no in edge weight) issue at root not seeing multiple components properly
    let -- preorder traversal
        middlePartList = concat $ fmap (getNewick fglGraph writeEdgeWeight writeNodeLable) (fmap (replicate 1) (G.out fglGraph index))
        label' = if writeNodeLable then label else T.empty -- trivial trees or write node name
    in
    -- trace ("MPL (" ++ (show $ length middlePartList) ++ ") " ++ show middlePartList ++ " " ++ show (G.out fglGraph index)) (
    -- "naked" root
    if null middlePartList then T.concat [T.singleton '(', label, T.singleton ')', T.singleton ';']
    -- single output edge
    else if length middlePartList == 1 then 
      T.concat [T.singleton '(', head middlePartList, T.singleton ')', label', T.singleton ';']
    else
      let middleText = T.intercalate (T.singleton ',') middlePartList
      in
      T.concat [T.singleton '(', middleText, T.singleton ')', label', T.singleton ';']
    --)


-- | makeLabel takes Maybe T.Text and retuns T.empty if Nothing, Text otherwise
makeLabel :: Maybe T.Text -> T.Text
makeLabel = fromMaybe T.empty

-- | fix for newick lack of paren in specific situation--inelegant
endStart :: T.Text
endStart = T.pack ")("

newEndStart :: T.Text
newEndStart = T.pack "),("

-- | getNewick takes an edge of a graph and either cretes the text if a leaf
-- or recurses down tree if has descendents, adding  commas, outper parens, labels, and edge weights if they exist.
-- need to filter redundant subtrees later at the forest level (Can have shared node between rooted components)
getNewick :: (Show a) => P.Gr T.Text a -> Bool -> Bool -> [G.LEdge a] -> [T.Text]
getNewick fglGraph writeEdgeWeight writeNodeLable inEdgeList
  | G.isEmpty fglGraph = [T.empty]
  | null inEdgeList = []
  | otherwise =
  let (_, curNodeIndex, edgeLabel) = head inEdgeList
      outEdges = G.out fglGraph curNodeIndex
  in
  -- is a leaf, no children
  if null outEdges then
    let leafLabel = G.lab fglGraph curNodeIndex
    in
    if isNothing leafLabel then error "Leaf without label in getNewick"
    else
      let newLabelList = if writeEdgeWeight then [T.concat [fromJust leafLabel, T.singleton ':', T.pack $ show edgeLabel]] else [fromJust leafLabel]
      in
      if length inEdgeList == 1 then newLabelList
      else [T.concat $ newLabelList ++ [T.singleton ','] ++ getNewick fglGraph writeEdgeWeight writeNodeLable (tail inEdgeList)]
  -- not a leaf, recurse
  else
    let nodeLabel = if not writeNodeLable then T.empty else makeLabel $ G.lab fglGraph curNodeIndex
        middlePartList = getNewick fglGraph writeEdgeWeight writeNodeLable (G.out fglGraph curNodeIndex)
    in
    if length middlePartList == 1 then  -- outdegree 1
      let middleText = T.replace endStart newEndStart (head middlePartList)
      in
      if not writeEdgeWeight then T.concat [T.singleton '(', middleText, T.singleton ')', nodeLabel]  : getNewick fglGraph writeEdgeWeight writeNodeLable  (tail inEdgeList)
      else T.concat [T.singleton '(', middleText, T.singleton ')', nodeLabel, T.singleton ':', T.pack $ show edgeLabel]  : getNewick fglGraph writeEdgeWeight writeNodeLable  (tail inEdgeList)
    else -- multiple children, outdegree > 1
      let middleText = T.intercalate (T.singleton ',') middlePartList
      in
      if not writeEdgeWeight then T.concat [T.singleton '(', middleText, T.singleton ')', nodeLabel]  : getNewick fglGraph writeEdgeWeight writeNodeLable  (tail inEdgeList)
      else T.concat [T.singleton '(', middleText, T.singleton ')', nodeLabel, T.singleton ':', T.pack $ show edgeLabel] : getNewick fglGraph writeEdgeWeight writeNodeLable  (tail inEdgeList)
      

-- |  stringGraph2TextGraph take P.Gr String a and converts to P.Gr Text a
stringGraph2TextGraph :: P.Gr String Double -> P.Gr T.Text Double
stringGraph2TextGraph inStringGraph =
    let (indices, labels) = unzip $ G.labNodes inStringGraph
        edges = G.labEdges inStringGraph
        textLabels = fmap T.pack labels
        newNodes = zip indices textLabels
    in
    G.mkGraph newNodes edges

    -- |  textGraph2StringGraph take P.Gr String a and converts to P.Gr Text a
textGraph2StringGraph :: P.Gr T.Text b -> P.Gr String b
textGraph2StringGraph inTextGraph =
    let (indices, labels) = unzip $ G.labNodes inTextGraph
        edges = G.labEdges inTextGraph
        stringLabels = fmap T.unpack labels
        newNodes = zip indices stringLabels
    in
    G.mkGraph newNodes edges

{-
Fucntions to relabel Dot greaph to RawGraph format
-}

-- | findStrLabel checks Attributes (list f Attribute) from Graphvz to extract the String label of node
-- returns Maybe Text
findStrLabel :: Attributes -> Maybe T.Text
findStrLabel = getFirst . foldMap getStrLabel

-- | getStrLabel takes an Attribute and reurns Text if StrLabel found, mempty otherwise
getStrLabel :: Attribute -> First T.Text
getStrLabel (Label (StrLabel txt)) = First . Just $ txt
getStrLabel _                      = mempty

-- | getLeafText takes a pairs (node vertex number, graphViz Attributes)
-- and returns Text name of leaf of Stringified nude number if unlabbeled
getLeafText :: (Int, Attributes) -> T.Text
getLeafText (nodeIndex, nodeLabel) =
  let maybeTextLabel = findStrLabel nodeLabel
  in
 fromMaybe (T.pack $ show nodeIndex)  maybeTextLabel

 -- | splitVertexList splits the vertices of a graph into ([root], [leaf], [non-leaf-non-root])
splitVertexList ::  P.Gr a b -> ([G.LNode a], [G.LNode a], [G.LNode a])
splitVertexList inGraph =
  if G.isEmpty inGraph then ([],[],[])
  else
    let -- leaves
        degOutList = G.outdeg inGraph <$> G.nodes inGraph
        newNodePair = zip degOutList (G.labNodes inGraph)
        leafPairList = filter ((==0).fst ) newNodePair
        (_, leafList) = unzip leafPairList

        -- roots
        degInList = G.indeg inGraph <$> G.nodes inGraph
        newRootPair = zip degInList (G.labNodes inGraph)
        rootPairList = filter ((==0).fst ) newRootPair
        (_, rootList) = unzip rootPairList

        -- non-leaves, non-root
        nodeTripleList = zip3 degOutList degInList (G.labNodes inGraph)
        nonLeafTripleList = filter ((>0).fst3 ) $ filter ((>0).snd3 ) nodeTripleList
        (_, _, nonLeafList) = unzip3 nonLeafTripleList   
    in
    (rootList, leafList, nonLeafList)

-- | getLeafList returns leaf complement of graph from DOT file
getLeafList ::  P.Gr Attributes Attributes -> [G.LNode T.Text]
getLeafList inGraph =
  if G.isEmpty inGraph then []
  else
    let degOutList = G.outdeg inGraph <$> G.nodes inGraph
        newNodePair = zip degOutList (G.labNodes inGraph)
        leafPairList = filter ((==0).fst ) newNodePair
        (_, leafList) = unzip leafPairList
        (nodeVerts, _) = unzip leafList
        newLabels = fmap getLeafText leafList
        leafList' = zip nodeVerts newLabels
    in
    leafList'

--  | relabelFGL takes P.Gr Attributes Attributes and converts to P.Gr T.Text Double
relabelFGL :: P.Gr Attributes Attributes -> P.Gr T.Text Double
relabelFGL inGraph =
  if G.isEmpty inGraph then G.empty 
  else 
    let newLeafList = getLeafList inGraph
        newEdgeLIst = fmap relabeLEdge (G.labEdges inGraph)
    in
    G.mkGraph newLeafList newEdgeLIst

-- | relabeLEdge convertes edhe labels to Double
relabeLEdge :: G.LEdge b -> G.LEdge Double
relabeLEdge (u,v,_) = (u,v,0.0:: Double)

--  | relabelFGL takes P.Gr Attributes Attributes and converts to P.Gr T.Text Double
relabelFGLEdgesDouble :: P.Gr a b -> P.Gr a Double
relabelFGLEdgesDouble inGraph =
  if G.isEmpty inGraph then G.empty 
  else  
    let newEdgeList = fmap relabeLEdge (G.labEdges inGraph)
    in
    G.mkGraph (G.labNodes inGraph) newEdgeList

-- | convertGraphToStrictText take a graphs with laze Text and makes it strict.
convertGraphToStrictText :: P.Gr T.Text Double -> P.Gr StrictT.Text Double
convertGraphToStrictText inGraph =
  if G.isEmpty inGraph then G.empty
  else 
    let nodeList = G.labNodes inGraph
        nodesStrictText = fmap T.toStrict $ fmap snd nodeList
        nodeIndices = fmap fst nodeList
    in
    G.mkGraph (zip nodeIndices nodesStrictText) (G.labEdges inGraph)


