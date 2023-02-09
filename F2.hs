import Data.List

----1-----
data MolSeq = DNA { name :: String , sequence :: String } | Protein { name :: String , sequence :: String } deriving (Show)  --"deriving Show" för att kunna printa


----2-----
string2seq :: String -> String -> MolSeq
string2seq name seq = if any (`elem` "RNDEQHILKMFPSWYV") seq
    then (Protein name seq)
    else (DNA name seq)

----3-----
--seqName :: MolSeq -> String
--seqName (_ a _) = a

seqSequence :: MolSeq -> String
seqSequence (DNA _ b) = b

seqLength :: MolSeq -> Int
seqLength (DNA _ b) = length b

----4-----
{-
seqDistance :: MolSeq -> MolSeq -> Double
seqDistance (DNA _ a) (DNA _ b) = -3.0/4.0*log(1.0-4.0*(findDiff(a,b)/fromIntegral(length(a)))/3.0)
seqDistance (Protein _ a) (Protein _ b) = -19.0/20.0*log(1.0-20.0*(findDiff(a,b)/fromIntegral(length(a)))/19.0)
seqDistance (DNA _ _) (Protein _ _) = error "Cannot compare DNA to Protein"
seqDistance (Protein _ _) (DNA _ _) = error "Cannot compare Protein to DNA"
-}


findDiff :: (String, String) -> Int
findDiff ([],[]) = 0
findDiff (x:xs, y:ys) = if (y==x)
    then findDiff (xs, ys)
    else 1 + findDiff (xs, ys)

-------3.1-----------
data Profile = Profile {matrix :: [[(Char,Double)]], profType :: String , numSeq :: Int, profName :: String} deriving (Show)

molseqs2profile :: String -> [MolSeq] -> Profile
molseqs2profile s ms = Profile m pt ns pn 
    where 
    m = map (map (\x -> (fst x, fromIntegral(snd x) / fromIntegral(length(ms))))) (makeProfileMatrix ms)
    pt = seqType(head ms) 
    ns = length(ms) 
    pn = s 

profileName :: Profile -> String
profileName (Profile _ _ _ pn) = pn

profileFrequency :: Profile -> Int -> Char -> Double
profileFrequency (Profile m _ _ _) i c = snd((findLetter (transpose m) c) !! i)

findLetter :: [[(Char, Double)]] -> Char -> [(Char, Double)]
findLetter [] c = []
findLetter (x:xs) a = if (fst(x!!0) == a)
    then x
    else findLetter xs a



seqType :: MolSeq -> String
seqType (DNA _ _) = "DNA"
seqType (Protein _ _) = "Protein"

------hjälpkod-------
nucleotides = "ACGT"
--aminoacids = sort "ARNDCEQGHILKMFPSTWYV"
aminoacids = "ACDEFGHIKLMNPQRSTVWY" 
makeProfileMatrix :: [MolSeq] -> [[(Char,Int)]] -- ändrat returtyp
makeProfileMatrix [] = error "Empty sequence list"
makeProfileMatrix sl = res
    where
        t = seqType (head sl)
        defaults =
            if (t == "DNA") then
                zip nucleotides (replicate (length nucleotides) 0) -- Rad (i)
            else
                zip aminoacids (replicate (length aminoacids) 0) -- Rad (ii)
        strs = map seqSequence sl -- Rad (iii)
        tmp1 = map (map (\x -> ((head x), (length x))) . group . sort) (transpose strs) -- Rad (iv)
        equalFst a b = (fst a) == (fst b)
        res = map sort (map (\l -> unionBy equalFst l defaults) tmp1)
    
    

main::IO()
main = do
    let mol1 = string2seq "mol1" "ACATAA"
    let mol2 = string2seq "mol2" "AAGTCA"
    let mol3 = string2seq "mol3" "ACGTGC"
    let mol4 = string2seq "mol4" "AAGTTC"
    let mol5 = string2seq "mol5" "ACGTAA"
    let profiltest = molseqs2profile "profiltest" [mol1, mol2, mol3, mol4, mol5]
    print(profiltest)
    print(profileName(profiltest))
    print(profileFrequency profiltest 4 'T')
