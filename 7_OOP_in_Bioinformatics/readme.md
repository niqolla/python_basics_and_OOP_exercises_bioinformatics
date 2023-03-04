1) Define a new class named Protein, with the following definition. If
necessary, you can define the private methods you need

        # Protein

        +identifier: String

        +sequence: String

        ---

        +get_identifier(): string

        +get_sequence(): string

        +get_mw(): float

        +has_subsequence( Protein): boolean

        +get_length(): integer


2) Modify the FASTA_iterator generator function to yield Protein objects
instead of tuples. In each iteration, the function must yield a Protein
Object:
FASTA_iterator( fasta_filename )