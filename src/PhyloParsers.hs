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

{-- to do
    Add fgl <-> Dot
    Add in Fastc
    TNT later (or simple for now) 
--}

{- 

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

module PhyloParsers (forestEnhancedNewickList2FGLList,
                     fglList2ForestEnhancedNewickString
                    ) where

import Data.Maybe
import qualified Data.Graph.Inductive.Graph as G
import qualified Data.Graph.Inductive.PatriciaTree as P
import qualified Data.Text.Lazy as T
import Data.Char (isSpace)
import qualified Data.Map.Strict as Map

{--  
    Using Text as ouput for non-standard ascii charcaters (accents, umlautes etc)


    ToDo:
      Parallelize
--}

-- | function for first element of triple 
fst3 :: (a,b,c) -> a
fst3 (d,_,_) = d

-- | getForestEnhancedNewickList takes String file contents and returns a list 
-- of fgl graphs with Text labels for nodes and edges or error if not ForestEnhancedNewick or Newick formats.
forestEnhancedNewickList2FGLList :: String -> [P.Gr T.Text Double]
forestEnhancedNewickList2FGLList fileString = 
    if null fileString then error "Empty file string input in getForestEnhancedNewickList"
    else 
        let fileText = T.pack fileString
            feNewickList = fmap removeNewickSpaces $ fmap removeNewickComments $ divideGraphText fileText
        in
        fmap text2FGLGraph feNewickList
 
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
            let firstPart = T.snoc ((T.takeWhile (/= ';')) inText) ';'
                restPart = T.tail $ (T.dropWhile (/= ';')) inText
            in
            firstPart : divideGraphText restPart
        else error "First character in graph representation is not either < or ("

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
  if T.null inText then inText
  else 
    let firstPart = T.replace (T.singleton ' ') (T.singleton '_') $ T.takeWhile (/= ''') (T.tail inText)
        restPart = T.tail $ T.dropWhile (/= ''') (T.tail inText)
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
    if firstChar == ''' then 
      let (newText, restTest) = convertQuotedText inText
      in
      T.concat [newText, removeNewickSpaces restText]  
    else if T.isSpace firstChar then removeNewickSpaces T.tail inText
    else T.cons firstChar (removeNewickSpaces T.tail inText)

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
        else if firstChar == '(' && lastChar == ';' then makeGraphFromPairList $ eNewick2FGL [] [] (-1, T.empty) [inGraphText] 
        else error ("Graph text not in ForestEnhancedNewick or (Enhanced)Newick format")


-- | fENewick2FGL takes a Forest Extended Newick (Text) string and returns FGL graph
-- breaks up forest and parses seprate eNewicks then modifes for any
-- common network nodes in the sub-graphs
fENewick2FGL :: T.Text -> P.Gr T.Text Double
fENewick2FGL inText =
  if T.null inText then error "Empty graph text in fENewick2FGL"
    else 
      -- split eNewicks
      let eNewickTextList = splitForest inText
      -- init to remove trailing ';' from eNewick
          eNewickGraphList = fmap makeGraphFromPairList $ fmap (eNewick2FGL [] [] (-1, T.empty) . (:[])) eNewickTextList
      in
      if length eNewickGraphList == 1 then head eNewickGraphList
      else
          -- merge graphs then merge network nodes and edges edges
          let fENewickGraph = mergeNetNodesAndEdges $ mergeFGLGraphs eNewickGraphList
          in
          fENewickGraph

-- | splitForest takes a Text (string) Forest Enhanced Newick representation and splits into 
-- its consituent Extended Newick representations
splitForest :: T.Text -> [T.Text]
splitForest inText = 
  if T.null inText then []
  else if (T.head inText /= '<') || (T.last inText /= '>') then error ("Invalid Forest Extended Newick representation," ++  
      " must begin with \'<\'' and end with \'>\' : " ++ (T.unpack inText))
  else 
    let partsList = filter (not.(T.null)) $ T.splitOn (T.singleton ';') inText
        eNewickList = fmap (T.append (T.singleton ';')) partsList
    in
    eNewickList

-- | makeGraphFromPairList takes pair of node list and edge list and returns Graph
-- | filters to remove place holder node and edges creted during eNewick pass
makeGraphFromPairList :: [(G.LNode T.Text,G.LEdge Double)] -> P.Gr T.Text Double
makeGraphFromPairList pairList = 
  if null pairLIst then G.empty
  else 
    let (nodeList, edgeList) = unZip pairList
    in
    G.mkGraph (filter ((> (-1)).fst) nodeList) (filter ((> (-1)).fst3) edgeList)

-- | getBranchLength extracts branch length from Text label and puts in '1' if there is no
-- branch length--makes sure after last ')'
getBranchLength :: T.Text -> Double
getBranchLength inText = 
  if (T.null inText) then error "Null text in getBranchLength"
  else 
    let a =  T.reverse $ T.takeWhile (/= ':') $ T.takeWhile (/= ')') T.reverse inText
    in
    if T.null a then 1
    else (read (T.unpack a) :: Double)

-- | getNodeLabel get--or makes--a label for a node
-- after last ')' before any ':', without ',' after last ')'
getNodeLabel :: Int -> T.Text -> T.Text
getNodeLabel nodeNumber inText =
  if T.null inText then error "Null text in getNodeLabel" 
  else 
    let a = T.takeWhile (/= ':') $ T.reverse $ T.takeWhile (/= ')') $ T.reverse inText
    in
    if (T.any (==',') a) then (T.append (T.pack $ show nodeNumber) (T.pack "HTU")) 
    else if (T.null a) then (T.append (T.pack $ show nodeNumber) (T.pack "HTU")) 
    else a
            
-- | getLeafInfo takes Text of teminal (no ',') and parses to yeild
-- either a single leaf label, edge, and edge weight, or two
-- leaves with labels and costs if there is a network node as parent
-- also searches for existing node labels to make correct edge and node 
-- if no new node is needed to be crearted (since already exists due to
-- network set dummy node with index -1 to be finltered later
-- scans for existing nodes and retuns if found, and sets appropriate edge indices
getLeafInfo :: T.Text -> G.LNode T.Text -> [G.LNode T.Text] -> [(G.LNode T.Text,G.LEdge Double)]
getLeafInfo leafText parentNode nodeList = 
  if T.null leafText then error "Empty leaf text in getLeafInfo"
  else 
    -- simple leaf
    if not (T.any (=='(') inText) then
      let leafLabel = T,takeWhile (/= ':') leafText
          edgeWeight = getBranchLength leafText
          -- CHECK FOR EXISTING
          thisNode = (length nodeList, leafLabel)
          thisEdge = (fst parentNode, length nodeList, edgeWeight)
          preexistingNode = checkForExistingNode leafLabel nodeList
      in
      if (preexistingNode == Nothing) then [(thisNode, thisEdge)]
      else 
        let newNode = fromJust preexistingNode
            newEdge = (fst parentNode, fst newNode, edgeWeight)
        in
        -- (-1) node index  filtered later to keep node list in good order
        [((-1, snd newNode), newEdge)]

    --complex leaf label
    --leaf and leafParent nodes
    -- sublabelText:Double
    else 
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
          leafNode = (1 + (length nodeList), leafLabel)
          leafEdge = (fst leafParentNode, fst leafNode, leafEdgeWeight)

          -- Check for existing nodes--assumes parent and child do not have same label
          existingLeafNode = checkForExistingNode leafLabel nodeList
          existingleafParentNode = checkForExistingNode leafParentLabel nodeList
      in
      if ((existingLeafNode == Nothing) && (existingleafParentNode = Nothing)) then [(leafNode, leafEdge),(leafParentNode, leafParentEdge)]
      else if (existingLeafNode == Nothing) then 
        -- leaf existing, parent new
        let newLeaf = fromJust existingLeafNode
            newLeafNode = (-1, snd newLeaf)
            newLeafEdge = (fst leafParentNode, fst newLeaf, leafEdgeWeight)
        in
        [(newLeafNode, newLeafEdge),(leafParentNode, leafParentEdge)]


      else if (existingleafParentNode == Nothing) then 
        -- parent existing, leaf new
        let newParentNode = fromJust existingleafParentNode
            newLeafParentNode = (-1, snd newParentNode)
            newLeafParentEdge = (fst parentNode, fst newParentNode, leafParentEdgeWeight)
        in
        [(leafNode, leafEdge),(newLeafParentNode, newLeafParentEdge)]

      else 
        -- both leaf and parent existed (not actually allowed in "phylogenetic" graph but not excluded by format)
        let newLeaf = fromJust existingLeafNode
            newParentNode = fromJust existingleafParentNode
            newLeafNode = (-1, snd newLeaf)
            newLeafParentNode = (-1, snd newParentNode)
            newLeafEdge = (fst newParentNode, fst newLeaf, leafEdgeWeight)
            newLeafParentEdge = (fst parentNode, fst newParentNode, leafParentEdgeWeight)
        in
        [(newLeafNode, newLeafEdge),(newLeafParentNode, newLeafParentEdge)]


-- | getBodyParts takes a Text of a subTree and splits out the group description '(blah)', any node label
-- and any branch length
getBodyParts :: T.Text -> Int -> (T.Text, T.Text, Double)
getBodyParts inRep nodeNumber = 
  if T.null inRep then error "No group to parse in getBodyParts"
  else 
      let subGraphPart =  T.reverse $ T.dropWhile (/= ')') $ T.reverse inRep
          branchLength =  getBranchLength inRep
          subGraphLabel = getNodeLabel nodeNumber inRep
      in
      (subGraphPart, subGraphLabel, branchLength)

-- | getChildren splits a subGraph Text '(blah, blah)' by commas, removing outer parens
getChildren :: T.Text -> [T.Text]
getChildren inText = 
  if T.null inText then []
  else if (T.head inText /= '(') || (T.last inText /= ')') then error ("Invalid Extended Newick component," ++  
      " must begin with \'(\'' and end with \')\' : " ++ (T.unpack inText))
  else 
    let guts = filter (not.(T.null)) $ T.splitOn (T.singleton ',') $ T.init $ T.tail inText
    in
    guts

-- | checkForExistingNode takes a node label and checs the node list for the first
-- node with the same label and returns a Maybe node, else Nothing
checkForExistingNode :: T.Text -> [G.LNode T.Text] -> Maybe (G.LNode T.Text)
checkForExistingNode nodeLabel nodeList =
  if null nodeList then Nothing
  else 
    let matchList = filter (==nodeLabel.snd) nodeList
    in
    if null matchList then Nothing
    else Just $ head matchList

-- | eNewick2FGL takes a single Extended Newick (Text) string and returns FGL graph
-- allows arbitrary in and out degree except for root and leaves
eNewick2FGL :: [G.LNode T.Text] -> [G.LEdge Double] -> G.LNode T.Text -> [T.Text] -> [(G.LNode T.Text,G.LEdge Double)]
eNewick2FGL nodeList edgeList parentNode inTextList = 
    if T.null inTextList then []
    else 
      let inText = head inTextList
      in  
      -- see if initial call and check format
      if null nodeList && ((T.head inText /= '(') || (T.last inText /= ';'))  then error ("Invalid Extended Newick component," ++  
      " must begin with \'(\'' and end with \')\' : " ++ (T.unpack inText))
      -- not first call and/or format OK
      else 
        let inText = T.takeWhile (/= ';') $ T.head inTextList  -- remove trailing ';' if first (a bit wasteful--but intial check on format)
        in
        -- is a single leaf
        if not (T.any (==',') inText) then 
          -- parse label ala Gary Olsen formalization
          -- since could have reticulate label yeilding two edges and two nodes
          -- Cardona et al 2008  Extended Newick
          getLeafInfo inText parentNode nodeList
        else 
          -- is subtree
          let (subTree, nodeLabel, edgeWeight) = getBodyParts inText (length nodeList)
              thisNode = (length nodeList, nodeLabel)
              thisEdge = (fst parentNode, length nodeList, edgeWeight)
              childTextList = getChildren subTree

              --check for existing node
              existingNode = checkForExistingNode nodeLabel nodeList
          in
          if existingNode == Nothing then (thisNode, thisEdge) : eNewick2FGL (thisNode : nodeList) (thisEdge : edgeList) thisNode childTextList
          else 
            let newNode = fromJust existingNode
                newEdge = (fst parentNode, fst newNode, edgeWeight)
            in
          -- allows the filtering out redundant nodes (-1 index) later, keeps node list in good shape for index determination
          ((-1, snd newNode)), newEdge) : eNewick2FGL nodeList (newEdge : edgeList) thisNode childTextList

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
mergeFGLGraphs curGraph inGraphList =
  if null inGraphList then curGraph
  else if G.isEmpty curGraph then mergeFGLGraphs (head inGraphList) (tail inGraphList)
  else 
    let firstGraph = head inGraphList
        firstNodes = G.labNodes firstGraph
        firstEdges = G.labEdges firstGraph
        curNodes = G.labNodes curGraph
        curEdges = G.labEdges curGraph
        newNodes = fmap (reindexNode (length curNodes)) firstNodes
        newEdges = fmap (reindexEdge (length curNodes)) firstEdges
    in
    mergeFGLGraphs (G.mkGraph (curNodes ++ newNodes) (curEgdes ++ newEdges)) (tail inGraphList)

-- | getNodeIndexPair take a list of unique nodes and checks successive nodes and 
-- adds to unique list, also creating a full list of pairs of indicess for non-unique that
-- can be used as an index map for edges
getNodeIndexPair :: [G.LNode T.Text] -> [G.LNode T.Text] -> ([G.LNode T.Text], [(Int, Int)])
getNodeIndexPair uniqueList pairList nodeToCheckList =
  if null nodeToCheckList then (reverse uniqueList, reverse pairList)
  else 
    let firstNode@(index, label) = head nodeToCheckList
        matchingNode = checkForExistingNode label uniqueList
    in
    if matchingNode == Nothing then getNodeIndexPair (firstNode : uniqueList) ((index, index) : pairList) (tail nodeToCheckList)
    else 
      let existingNode = fromJust matchingNode
          newPair = (index, fst existingNode)
      in
      getNodeIndexPair (uniqueList) (newPair : pairList) (tail nodeToCheckList)

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

-- | mergeNetNodesAndEdges takes a single graph and merges 
-- nodes and edges due to network nodes and edges
-- uses checkForExistingNode and creates a map from nodes to reindex edges
-- needs to be merged first if graphs are combined--or indices will be wrong
mergeNetNodesAndEdges :: P.Gr T.Text Double -> P.Gr T.Text Double
mergeNetNodesAndEdges inGraph = 
  if G.isEmpty inGraph then G.empty
  else 
    let nodeList = G.labNodes inGraph
        (uniqueNodeList, nodeIndexPairs) = getNodeIndexPair [] [] nodeList
        nodeMap = Map.fromList nodeIndexPairs
        reIndexedEdgeList = fmap (reIndexLEdge nodeMap) (G.labEdges inGraph)
    in
    G.mkGraph uniqueNodeList reIndexedEdgeList

-- | subTreeSize takes a nodeList and retuns the number of leaves that can be
-- traced back to those nodes (for single just pass list of 1 node)
subTreeSize :: P.Gr a a -> [G.LNode] -> Int -> Int 
subTreeSize inGraph nodeList counter =
  if null nodeList then counter
  else 
    let firstNode = head nodeList
        children = G.lsuc inGraph index
    in
    subTreeSize inGraph (tail nodeList ++ children) (counter + (length children))

-- | getRoot takes a greaph and list of nodes and returns vertex with indegree 0
-- so assumes a connected graph--with a single root--not a forest
getRoots :: P.Gr a a -> [G.LNode] -> [G.LNode]
getRoots inGraph nodeList =
  if null nodeList then [] --error "Root vertex not found in getRoot"
  else
    let firstNode@(index, label) = head nodeList
    in
    if (G.indeg inGraph index == 0) && (G.outdeg inGraph index > 0) then firstNode : getRoots inGraph (tail nodeList)
    else getRoots inGraph (tail nodeList)


-- | fgl2FEN take a fgl graph and returns a Forest Enhanced Newick Text
fgl2FEN :: (Show a) => P.Gr T.Text a -> Bool -> T.Text
fgl2FEN fglGraph writeEdgeWeight =
  if G.isEmpty fglGraph then error "Empty graph to convert in fgl2FEN"
  else 
    -- get forest roots
    let numRoots = getRoots inGraph (G.labNodes fglGraph)
        fenTextList = fmap (component2Newick fglGraph writeEdgeWeight) (sortOn (subTreeSize 0) numRoots)
        wholeRep = T.concat $ fmap (T.append (T.singleton '\n')) $ fenTextList
    in
    if (length fenTextList == 1) then wholeRep -- just a single tree/network
    else (T.snoc '>' $ T.cons '<' wholeRep)


-- | fglList2ForestEnhancedNewickString takes FGL representation of forest and returns
-- list of Forest Enhanced Newick as a single String
fglList2ForestEnhancedNewickString :: (Show a) => [P.Gr T.Text a] -> Bool -> String 
fglList2ForestEnhancedNewickString inFGLList writeEdgeWeight =
  if null inFGLList then "\n"
  else 
    let forestTextList = fmap (T.append (T.singleton '\n')) $ fmap fgl2FEN inFGLList
        forestListString = T.pack $ T.concat forestTextList
    in
    forestListString

-- | component2Newick take a graph and root and creartes enhanced newick from that root
component2Newick :: P.Gr T.Text a -> Bool -> G.LNode -> T.Text
component2Newick fglGraph writeEdgeWeight rootNode =
  if G.isEmpty fglGraph then error "Empty graph to convert in component2Newick"
  else 
    -- preorder pass
    T.empty