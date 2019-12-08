#!/usr/bin/env python

"""
    pyamrfinder images.
    Created Nov 2019
    Copyright (C) Damien Farrell

    This program is free software; you can redistribute it and/or
    modify it under the terms of the GNU General Public License
    as published by the Free Software Foundation; either version 3
    of the License, or (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program; if not, write to the Free Software
    Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
"""

try:
     import tkinter as tk
except:
     import Tkinter as tk

def logo_image():
     img = tk.PhotoImage(format='gif',data=
             'R0lGODlhQAFAAfcAAAAAAAwCAAYHCAkKDAgJBxsAABQGAA0OEAcMFhIUFhYY'
            +'GxkbHxQTCicAADcAACUXATQdBDgpADQoBwoWJxwfIw8dNRAhOyAiJyYpLyQn'
            +'KywwNi8zOjI3Pjk6OiotMkgAAFkAAGcAAHgAAFAyAm06BEk4KFtGA3pcAW1P'
            +'Bn9gAH1hCHthEntpLWVUMBQnRxgxWDk+RjY6Qh06aD5DTSJEekRKVElPWUxS'
            +'XVVUU3ZtT09VYVRbZ1hgbV1lc2VtfGlscnV1dWpeUYwAAIMAAJQAAJkAAIdf'
            +'AZhaBIxLCKdcBbRfBrFfCYFhA4JkCoxhA5dhBYVqE4dtG4hvH5RqGKVgBbRi'
            +'C7FgBrVmE7ZqG6xnGIpyJI12Lo96NYl1NJF+PJV3MZFoMbZtIrdxKLhyLLh1'
            +'Mrh6PLR4O6xuK4ZZLLl+Qqx6SZV5W5WFTJOBRJuNXJeJVI6CVbqDTLqARruJ'
            +'WKaBRpyPYZ+VbJ2SZY+LdLyOY72SaqSce6GYc72Wcr+bfbKSd7KPcSZLhypU'
            +'ly5cpTNltmpzgm54iHJ7jH9+hXV/kDZtxDlxzTx42Dt11H6JnHqEln6EhoGM'
            +'n4yMjaWfgL+fg6ijiL+hia2tnKqokbWjmoiTp4uWq4qWqpGds5qdqZSht5ej'
            +'upqnvpiht5egqaanpq+xo7K1rLGyp7e7ura7uaqrtMCkj8GrncGnlMGuocOy'
            +'q8GwpcO3s8S7vJ2qwqCuxaSyyrm/w6u507C/2bC91sW/w7vEy7rCxb/J1bPC'
            +'3b3H0rbF4bnJ5bzM6MXDy8bDxsXH0sTM3MbI1tbW1sXO4sHP6MXQ5cPR68jU'
            +'7MXW8snW8MnY9Mna+Ojo6PPz8wAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA'
            +'AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA'
            +'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')

     return img