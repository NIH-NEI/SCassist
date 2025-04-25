## Measuring Sentence Similarity with BERT

This report documents the code that was used to evaluate the semantic similarity between the input enrichment terms and the LLM response content, as part of the SCassist package. The report delves into the fascinating world of semantic similarity, exploring how the powerful BERT (Bidirectional Encoder Representations from Transformers) model can be harnessed to quantify the meaning shared between sentences. We'll dissect the code, providing a comprehensive understanding of each step and its significance.

**Narrative style code Documentation by :** Gemini \
**Documentation reviewed by :** Vijay Nagarajan PhD

**1. Setting the Stage: Libraries and Environment**

Our journey begins with importing the necessary tools:

```python
from transformers import BertTokenizer, BertModel
import torch
from sklearn.metrics.pairwise import cosine_similarity
import random
```

Before we can embark on our analysis, we need to set up a suitable environment:

- **Conda:** A package manager, conda helps us create and manage isolated environments for our project.
- **Environment Creation:**
   ```bash
   conda create -n mybert python=3.9.15 transformers
   ```
- **Activation:**
   ```bash
   conda activate mybert
   ```

**2. The Power of BERT: Loading the Model and Tokenizer**

At the heart of our analysis lies BERT, a transformer-based model renowned for its ability to capture the nuances of language. We load the pre-trained BERT model and its accompanying tokenizer:

```python
# Set the random seed
torch.manual_seed(42)
random.seed(42)

# Load the BERT tokenizer and model
tokenizer = BertTokenizer.from_pretrained('bert-base-uncased')
model = BertModel.from_pretrained('bert-base-uncased')
```

**3. Our Subjects: Example Sentences**

To illustrate the process, we'll use a simple example:

```python
# Example sentences
# sentence1 = "I like coding in Python."
# sentence2 = "Python is my favorite programming language."

# Read content from text files
with open("sentence1.txt", "r") as f:
    sentence1 = f.read().strip()

with open("sentence2.txt", "r") as f:
    sentence2 = f.read().strip()
```

**4. Tokenization: Breaking Down the Sentences**

The tokenizer breaks down the sentences into individual tokens:

```python
# Tokenize the sentences
tokens1 = tokenizer.tokenize(sentence1)
tokens2 = tokenizer.tokenize(sentence2)
```

**5. Special Tokens: Guiding the Model**

We add special tokens to the beginning and end of each sentence:

```python
# Add [CLS] and [SEP] tokens
tokens1 = ['[CLS]'] + tokens1 + ['[SEP]']
tokens2 = ['[CLS]'] + tokens2 + ['[SEP]']
```

**6. Chunking: Handling Long Sentences**

BERT has a maximum input length, so we split long sentences into chunks:

```python
# Segment the tokens into chunks
max_length = 512
chunk_size = max_length - 2  # Account for [CLS] and [SEP]
chunks1 = [tokens1[i:i + chunk_size] for i in range(0, len(tokens1), chunk_size)]
chunks2 = [tokens2[i:i + chunk_size] for i in range(0, len(tokens2), chunk_size)]
```

**7. Processing Chunks: Extracting Meaning**

We process each chunk through the BERT model, extracting its meaning:

```python
# Process each chunk
embeddings1 = []
embeddings2 = []
for chunk1, chunk2 in zip(chunks1, chunks2):
    # Add [CLS] and [SEP] tokens
    chunk1 = ['[CLS]'] + chunk1 + ['[SEP]']
    chunk2 = ['[CLS]'] + chunk2 + ['[SEP]']

    # Convert tokens to input IDs
    input_ids1 = torch.tensor(tokenizer.convert_tokens_to_ids(chunk1)).unsqueeze(0)
    input_ids2 = torch.tensor(tokenizer.convert_tokens_to_ids(chunk2)).unsqueeze(0)

    # Obtain the BERT embeddings
    with torch.no_grad():
        outputs1 = model(input_ids1)
        outputs2 = model(input_ids2)
        embeddings1.append(outputs1.last_hidden_state[:, 0, :])
        embeddings2.append(outputs2.last_hidden_state[:, 0, :])
```

**8. Concatenation: Combining the Embeddings**

We combine the embeddings from all chunks of each sentence:

```python
# Concatenate the embeddings from all chunks
embeddings1 = torch.cat(embeddings1, dim=0)
embeddings2 = torch.cat(embeddings2, dim=0)
```

**9. Cosine Similarity: Measuring the Angle**

We calculate the cosine similarity between the concatenated embeddings:

```python
# Calculate similarity
similarity_score = cosine_similarity(embeddings1, embeddings2)
print("Similarity Score:", similarity_score)
```

**Conclusion**

This report has illuminated the process of measuring semantic similarity using BERT. By tokenizing sentences, adding special tokens, chunking long sentences, and extracting embeddings, we can quantify the meaning shared between two sentences. This powerful technique opens doors to a wide range of applications, from text summarization and question answering to sentiment analysis and document classification. As we continue to explore the depths of natural language processing, BERT stands as a beacon of progress, guiding us towards a deeper understanding of the intricate tapestry of human language.
