import numpy as np
import sounddevice as sd


def generate_tone(char, duration_ms=500, sampling_freq=44100):
    """Generate a tone for a given character"""

    # Define a dictionary that maps characters to frequencies (in Hertz)
    # A frequence of 0 will be returned for an empty space character
    char_to_frequency = {
        'A': 440, 'B': 494, 'C': 523, 'D': 587, 'E': 659,
        'F': 698, 'G': 784, 'H': 880, 'I': 988, 'J': 1047,
        'K': 1175, 'L': 1319, 'M': 1397, 'N': 1568, 'O': 1760,
        'P': 1976, 'Q': 2093, 'R': 2349, 'S': 2637, 'T': 2794,
        'U': 3136, 'V': 3520, 'W': 3951, 'X': 4186, 'Y': 4699, 'Z': 5274,
        ' ': 0 
        }

    frequency = char_to_frequency.get(char.upper(), 0)
    if frequency > 0:
        n_samples = int(sampling_freq * duration_ms / 1000)
        t = np.linspace(0, duration_ms / 1000, n_samples, False)
        tone = 0.5 * np.sin(2 * np.pi * frequency * t)
        return tone
    else:
        print(f"Character '{char}' not found in the dictionary.")
        return None


if __name__ == "__main__":
    tone = generate_tone("Z")
    if tone is not None:
        sd.play(tone, samplerate=44100)
        sd.wait()
    else: 
        print("No tone generated")
