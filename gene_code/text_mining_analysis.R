rm(list = ls())

library(quanteda)
library(stringr)
library(R.matlab)

ab_file <- read.csv("M:/Desktop/Rcode/PPI_hub/S1_merge.csv")
abstracts <- ab_file$abstract

# ------ preprocessing and construct quanteda token ------
qcorpus <- corpus(abstracts)
tokens <- tokens(qcorpus, remove_punct = TRUE, remove_symbols = TRUE, remove_numbers = TRUE)
tokens <- tokens_tolower(tokens)
tokens <- tokens_remove(tokens, stopwords("english"))
tokens <- tokens_wordstem(tokens) 

# ------ Keyword dictionary (Supplementary Table 6) ------
keyword_dict <- dictionary(list(
  Verbal_communication = c("verbal*", "communicat* defic*", "communicat* impair*"),
  Language_speech = c("atypical language","delay* language","language comprehension",
                      "language* delay*","language* impair*","atypical speech",
                      "speech atypical*","speech development","speech impediment*",
                      "dysphonia","prosod*","voice control","vocal cord","vocal product*"),
  Social = c("social adjustment","social affect","social behavi*","social bonding",
             "social cognit*","social competence","social deficit*","social development",
             "social impair*","social interact*","social test","social withdraw*"),
  Attention_deficits = c("ADHD","attention*deficit","attention problem*","inattention"),
  Repetitive_restricted = c("hand*flapping","repetitive behavi*","stereoty* action*",
                            "stereoty* behav*","stereotypy","stereoty* activit*",
                            "stereotyp* mo*","stimming","circumscribed interest*",
                            "hyper*focus*","inflexib*","persev* interest*","restrict* behavi*",
                            "restrict* interest*","restricted*repetitive","ritualistic","sameness"),
  Compulsive_impulsive_obsessive = c("compulsiv*","*grooming","impulsiv*","motor impulsiv*",
                                     "motor tic*","obsess*"),
  Self_harm = c("completed suicide","self*harm*","self*injur*","suicidal behav*",
                "suicidal idea*","suicidal thought*","suicidality")
))

dfm <- dfm(tokens)
dfm_lookup <- dfm_lookup(dfm, dictionary = keyword_dict)
print(dfm_lookup)

# ------ calculate frequency ------
keyword_counts <- colSums(dfm_lookup)
relative_freq <- keyword_counts / sum(keyword_counts)

print("count：")
print(keyword_counts)
print("frequency：")
print(relative_freq)

writeMat("PLSR_result/output_gene_order_fromGLM/PPI_genes/text_mining/S1.mat",freq = relative_freq)
