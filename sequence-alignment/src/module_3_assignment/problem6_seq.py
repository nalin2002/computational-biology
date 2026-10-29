from pathlib import Path
from collections import defaultdict
from time import perf_counter


def clean_text(text):
    return "".join(char.lower() for char in text if char.isalpha())


def load_text(file_path):
    with open(file_path, "r", encoding="utf-8") as file:
        return clean_text(file.read())


def load_words(file_path):
    with open(file_path, "r", encoding="utf-8") as file:
        return {
            word
            for line in file
            if len(word := clean_text(line.strip())) >= 4
        }


def build_word_index(all_words):
    word_index = defaultdict(list)

    for word in all_words:
        word_index[word[0]].append(word)

    return word_index


def find_words(els_text, word_index, names):
    for pos, char in enumerate(els_text):
        for word in word_index.get(char, []):
            if els_text.startswith(word, pos):
                yield pos, word, word in names


def find_phrases(entity_matches, window=25):
    phrases = []
    cluster = []

    for pos, word in entity_matches:
        if not cluster:
            cluster = [(pos, word)]

        elif pos - cluster[0][0] <= window:
            cluster.append((pos, word))

        else:
            if len(cluster) >= 2:
                phrases.append(cluster)

            cluster = [(pos, word)]

    if len(cluster) >= 2:
        phrases.append(cluster)

    return phrases


def process_book(text, book_name, word_index, names):
    output_dir = Path("../../output/biblecode/benchmark")
    output_dir.mkdir(parents=True, exist_ok=True)

    words_path = output_dir / f"{book_name}-words.txt"
    phrases_path = output_dir / f"{book_name}-phrases.txt"

    with open(words_path, "w", encoding="utf-8") as words_file, \
         open(phrases_path, "w", encoding="utf-8") as phrases_file:

        # Strides 2 to 100 inclusive
        for stride in range(2, 101):
            els_text = text[::stride]
            entity_matches = []

            for pos, word, is_entity in find_words(
                els_text, word_index, names
            ):
                words_file.write(f"{stride},{pos},{word}\n")

                if is_entity:
                    entity_matches.append((pos, word))

            phrases = find_phrases(entity_matches, window=25)

            for phrase in phrases:
                values = [str(stride)]

                for pos, word in phrase:
                    values.extend([str(pos), word])

                phrases_file.write(",".join(values) + "\n")


# Load reference words
words = load_words("../../data/words_10k.txt")
names = load_words("../../data/names_and_dates.txt")

all_words = words | names
word_index = build_word_index(all_words)

# Load Moby Dick
text = load_text("../../data/moby-dick.txt")

# Measure sequential execution time
start = perf_counter()

process_book(text, "moby-dick", word_index, names)

end = perf_counter()

print(f"Moby Dick (strides 2–100): {end - start:.4f} seconds")